#include <cstring>
#include "raylist.h"
#include "communicator.h"
#include "buffer.h"

Communicator::Communicator() {
    MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    master = size - 1;
    first = true;
    
    int group = rank == master ? 1 : 0;
    MPI_Comm_split(MPI_COMM_WORLD, group, rank, &Client_Comm);
    MPI_Comm_size(Client_Comm, &group_size);
    MPI_Comm_rank(Client_Comm, &group_rank);
    
    compress_buffer.reserve(1048608 * 21);
    quit = false;
    pause = false;

    os = std::ofstream("out/proc_" + std::to_string(rank));
}
    
Communicator::~Communicator() {
    MPI_Finalize();
}

void Communicator::Isend_rays(struct Rays* buffer, int size, int dst, int tag) {
    int buffer_size = buffer->get_size();
    int width = buffer->store_width;
    float *data = buffer->get_data();
    int    send_size = buffer_size > size ? size : buffer_size;
 //   printf("send ray send size %d buffer_size%d size %d\n", send_size, buffer_size, size);    
    MPI_Isend(&data[(buffer_size - send_size) * width], send_size * width, MPI_FLOAT, dst, 1, MPI_COMM_WORLD, &req[tag]);
    buffer->size -= send_size;
}

void Communicator::mpi_wait(int tag){
     MPI_Wait(&req[tag],sta);
}

void Communicator::send_noray(int dst) {
//    printf("send norays %d %d\n", rank, dst);
    int msg[MSG_SIZE];
    msg[0] = 0;
    MPI_Send(&msg[0], MSG_SIZE, MPI_INT, dst, 1, MPI_COMM_WORLD);
}

void Communicator::send_end(int dst) {
    int msg[MSG_SIZE];
    msg[0] = -1;
    MPI_Send(&msg[0], MSG_SIZE, MPI_INT, dst, 1, MPI_COMM_WORLD);
}

void Communicator::reduce_image(float* film, float *reduce_buffer, int pixel_num, bool server){
    if(server){
        MPI_Reduce(film, reduce_buffer, pixel_num, MPI_FLOAT, MPI_SUM, 0, Client_Comm); 
    } else {
        MPI_Reduce(film, reduce_buffer, pixel_num, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
    }
}

void Communicator::send_msg(int dst, int* msg){
    MPI_Send(&msg[0], MSG_SIZE, MPI_INT, dst, 1, MPI_COMM_WORLD);
}

void Communicator::recv_msg(int dst, int *msg) {
    MPI_Recv(msg, MSG_SIZE, MPI_INT, dst, 1, MPI_COMM_WORLD, sta);
}

void Communicator::recv_rays(int src, int recv_size, struct Rays* raylist) {
    if(recv_size == 0) return;
    int width = raylist->store_width;
    os<< rank << " <- " << src << "|recv" <<  recv_size <<" "<<width<<std::endl;
    float *rays = &raylist->get_data()[raylist->get_size() * width];
    raylist->size += recv_size;

    MPI_Recv(rays, recv_size * width, MPI_FLOAT, src, 1, MPI_COMM_WORLD, sta); 

//    for(int i = 0; i < 10; i ++) {
//        printf("|r %d %d %d", ids[i * width], ids[i * width + 9], ids[i * width + 15]);
//    }
//    printf("\n");
}

void Communicator::send_rays(int dst, int send_size, struct Rays* buffer) {
    if(send_size == 0) return;
    int width = buffer->store_width;
    float *rays = &buffer->get_data()[(buffer->get_size() - send_size) * width];
    buffer->size -= send_size;
    os << rank << "-> " <<  dst << "|send " << send_size << width << std::endl;
    MPI_Send(rays, send_size * width, MPI_FLOAT, dst, 1, MPI_COMM_WORLD);
}



void Communicator::purge_completed_mpi_buffers() {
    bool done = false;
    while (! done) {
        done = true;
        for (std::vector<mpi_send_buffer *>::iterator i = mpi_in_flight.begin(); done && i != mpi_in_flight.end(); i++)
        {
            int lflag, rflag; MPI_Status s;
            mpi_send_buffer *m = (*i);
      
            MPI_Test(&m->lrq, &lflag, &s);
      
            if (m->n > 1)
                MPI_Test(&m->rrq, &rflag, &s);
            else
                rflag = true;
               
            if (lflag && rflag) // if left is gone or, if bcast, BOTH are gone, then
            {
                done = false;
       
                mpi_in_flight.erase(i);
                free(m->send_buffer);
       
                delete m;
            }
        }
    }
}

int Communicator::Export(Message *m, ProcStatus *rs) {
	int root = m->get_root();
	static int t = 0;
	int k = 0;
	// My rank relative to the broadcast root
	int d = ((size + rank) - root) % size;
    
    os << "mthread root " << root<< " rank " << rank << std::endl;
	// Only export if its either P2P or broadcast and this isn't a leaf
	if (m->is_broadcast() && (2*d + 1) >= size)
		return 0;

	int tag = (MPI_TAG_UB) ? t % MPI_TAG_UB : t % 65535;
	t++;

	// If its a broadcast message, choose up to two destinations based
	// on the broadcast root, the rank and the size.  Otherwise, just ship it.
    
	struct mpi_send_buffer *msb = new mpi_send_buffer;
    
    msb->total_size = m->get_header_size() + (m->has_content() ? m->get_size() : 0);
    msb->send_buffer = (char *)malloc(msb->total_size);
    memcpy(msb->send_buffer, m->get_header(), m->get_header_size());
    os<< "mthread msb->total size " << msb->total_size << " " <<m->get_size()<<"broadcast"<<m->is_broadcast()<<std::endl;
    if (m->has_content())
        memcpy(msb->send_buffer + m->get_header_size(), m->get_content(), m->get_size());
    
    if (m->is_broadcast()) {
        os<<"mthread broadcast\n";
        int l = (2 * d) + 1;
		int destination = (root + l) % size;
		MPI_Isend(msb->send_buffer, msb->total_size, MPI_CHAR, destination, tag, MPI_COMM_WORLD, &msb->lrq);
		k++;

		if ((l + 1) < size)
		{
			destination = (root + l + 1) % size;
		    MPI_Isend(msb->send_buffer, msb->total_size, MPI_CHAR, destination, tag, MPI_COMM_WORLD, &msb->rrq);
			k++;
		}
    } else {
        os<< "mthread send rays" << m->get_ray_size() <<"to "<<m->get_destination()<< std::endl;
        MPI_Isend(msb->send_buffer, msb->total_size, MPI_CHAR, m->get_destination(), tag, MPI_COMM_WORLD, &msb->lrq);
		k++;
	}
    if(m->has_content())
        rs->accumulate_sent(k * m->get_ray_size());
	msb->n = k;
	mpi_in_flight.push_back(msb);
	
    return k;
}

void Communicator::collective(ProcStatus *rs) {

    rs->lock();

    int global_counts[3]; // [4] -> [3] camera active is not necessary
	int local_counts[] = {rs->all_thread_waiting() ? 0 : 1, rs->get_sent_ray_count(), rs->get_recv_ray_count()};

    MPI_Allreduce(local_counts, global_counts, 3, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    os << "mthread all end collective: local " << local_counts[0] << "sent " << local_counts[1] <<"recv "<< local_counts[2]<< std::endl;
    os << "mthread all end collective: global " << global_counts[0] << "sent " << global_counts[1] <<"recv "<< global_counts[2]<< std::endl;
	// If no raylists exist anywhere, and the number of received pixels matches the number of sent pixels,
	// and there are no camera rays currently being generated, this rendering is done.

//    rs->unlock();
//    return true;
     
    if (global_counts[0] == 0 && global_counts[1] <= 10000 + global_counts[2]) {
        //all process done, no local rays, global sent rays == recv rays 
        rs->set_exit();
        rs->unlock();
    } else {
        if (global_counts[0] != 0)
             std::cerr << "thread havent finish"<<global_counts[0]<<" "<<global_counts[1]<<"\n";
        else if (global_counts[1] != global_counts[2])
             std::cerr << "rays havent recv"<<global_counts[0]<<"\n";
        os<<"set proc busy\n"; 
        rs->set_proc_busy(rank);
        rs->unlock();
    }
}

bool Communicator::recv_message(RayList** List, ProcStatus *rs) {
    int recv_ready;
    MPI_Status status;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_ready, &status);

    if (recv_ready) {
        Message *incoming_message = new Message(status, MPI_COMM_WORLD); 
        os<<"mthread recv rays from "<<incoming_message->get_sender()<<"size "<<incoming_message->get_ray_size()
          <<"root: "<<incoming_message->get_root()<<"chunk "<<incoming_message->get_chunk()<<"\n";
        if (incoming_message->is_broadcast()) {
            os << "mthread recv broadcast\n";
            Export(incoming_message, rs); 
            if (incoming_message->is_collective()) {
                //check if need to synchronous
                collective(rs);
            } else {
                os << "mthread| set proc " << incoming_message->get_root() << "idle \n";
                rs->set_proc_idle(incoming_message->get_root());
            }
        } else {
            os << "mthread| incoming chunk "<<incoming_message->get_chunk()<< "root" << incoming_message->get_root() << "\n";
            RayList *in = List[incoming_message->get_chunk()];
            rs->accumulate_recv(incoming_message->get_ray_size());
            incoming_message->deserialize(in);
            os << "mthread| after deserialize "<<in->size()<<std::endl;
            //set itself busy
            rs->set_proc_busy(rank);
        }
        return true; 
    } else {
        return false;
    } 
}

bool Communicator::send_message(Message *outgoing_message, ProcStatus *rs) {

    //outList lock
    //serialize
    
    //current process is not destination or broadcast, send to destination 
    os << "mthread| send message"<<outgoing_message->get_destination() <<"\n";
    if (outgoing_message->is_broadcast() || (outgoing_message->get_destination() != rank))
        Export(outgoing_message, rs);
    if (outgoing_message->is_broadcast()) {
        if (outgoing_message->is_collective()) {
            os << "mthread|send collective\n";    
            collective(rs);
            // We copy the message in Export, and put the copy on the list to be held until 
            // the message actually leaves.   So here we acknowlege that the
            // collective action finishes.   The collective action can block, so we won't get here
            // until the messages have arrived down the tree

//            if (outgoing_message->isBlocking())
//            {
//              outgoing_message->Signal();        // blocked guy will delete
//            }
        } else {
            os << "mthread| broadcast status";    
        } 
//        else {
//            // If we enqueue the message for async work, we will wait until its been performed
//            // (in MessageManager::workThread) to signal or delete the message
//            
//            income_queue->enqueue(outgoing_message);
//            message_deserialize(inList);
    }
    
    return quit;
}

bool Communicator::barrier(struct RayList* out) {

}

