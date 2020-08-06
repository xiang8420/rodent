#include <cstring>
#include "RayList.h"
#include "communicator.h"
#include "../driver/buffer.h"

Communicator::Communicator() {
    MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    master = size - 1;
    
    pure_master = (size > 1 && size % 2 == 1); // if master join rendering 
    
    int group = rank == master ? 1 : 0;
    MPI_Comm_split(MPI_COMM_WORLD, group, rank, &Client_Comm);
    MPI_Comm_size(Client_Comm, &group_size);
    MPI_Comm_rank(Client_Comm, &group_rank);
    
    printf("comm rank %d \n", rank);
    os = std::ofstream("out/proc_" + std::to_string(rank));
    send_ray_count = 0; recv_ray_count = 0;
    send_msg_count = 0; recv_msg_count = 0;
}
    
Communicator::~Communicator() {
    printf("delete communicator\n");
    std::cout<<"msg send "<< send_msg_count<<"recv "<<recv_msg_count
      <<"ray send "<< send_ray_count<<"recv "<<recv_ray_count<<"\n";
    MPI_Finalize();
}

void Communicator::Isend_rays(struct Rays* buffer, int size, int dst, int tag) {
    int buffer_size = buffer->get_size();
    int width = buffer->store_width;
    float *data = buffer->get_data();
    int send_size = buffer_size > size ? size : buffer_size;
 //   printf("send ray send size %d buffer_size%d size %d\n", send_size, buffer_size, size);    
    MPI_Isend(&data[(buffer_size - send_size) * width], send_size * width, MPI_FLOAT, dst, 1, MPI_COMM_WORLD, &req[tag]);
    buffer->size -= send_size;
    send_ray_count ++;
}

void Communicator::mpi_wait(int tag){
     MPI_Wait(&req[tag],sta);
}

void Communicator::send_noray(int dst) {
//    printf("send norays %d %d\n", rank, dst);
    int msg[MSG_SIZE];
    msg[0] = 0;
    MPI_Send(&msg[0], MSG_SIZE, MPI_INT, dst, 1, MPI_COMM_WORLD);
    send_msg_count++;
}

void Communicator::send_end(int dst) {
    int msg[MSG_SIZE];
    msg[0] = -1;
    MPI_Send(&msg[0], MSG_SIZE, MPI_INT, dst, 1, MPI_COMM_WORLD);
    send_msg_count++;
}

void Communicator::reduce_image(float* film, float *reduce_buffer, int pixel_num){
    if(pure_master){
        MPI_Reduce(film, reduce_buffer, pixel_num, MPI_FLOAT, MPI_SUM, 0, Client_Comm); 
    } else {
        MPI_Reduce(film, reduce_buffer, pixel_num, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
    }
}

void Communicator::send_msg(int dst, int* msg){
    MPI_Send(&msg[0], MSG_SIZE, MPI_INT, dst, 1, MPI_COMM_WORLD);
    send_msg_count++;
}

void Communicator::recv_msg(int dst, int *msg) {
    MPI_Recv(msg, MSG_SIZE, MPI_INT, dst, 1, MPI_COMM_WORLD, sta);
    recv_msg_count++;
}

void Communicator::recv_rays(int src, int recv_size, struct Rays* raylist) {
    if(recv_size == 0) return;
    int width = raylist->store_width;
    os<< rank << " <- " << src << "|recv" <<  recv_size <<" "<<width<<std::endl;
    float *rays = &raylist->get_data()[raylist->get_size() * width];
    raylist->size += recv_size;

    MPI_Recv(rays, recv_size * width, MPI_FLOAT, src, 1, MPI_COMM_WORLD, sta); 

    recv_ray_count++;
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
    send_ray_count++;
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
    
    os << "mthread root " << root<< " rank " << rank <<" 2*d + 1 "<<2*d+1<< std::endl;
	// Only export if its either P2P or broadcast and this isn't a leaf
	if (m->is_broadcast() && (2*d + 1) >= size) {
        os << "mthread root | " << root<< " rank " << rank <<" 2*d + 1 "<<2*d+1<< std::endl;
        return 0;
    }

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
     
    os<<"mthread comm after copy content\n";
    if (m->is_broadcast()) {
        os<<"mthread broadcast\n";
        int l = (2 * d) + 1;
		int destination = (root + l) % size;
		MPI_Isend(msb->send_buffer, msb->total_size, MPI_CHAR, destination, tag, MPI_COMM_WORLD, &msb->lrq);
		k++;

		if ((l + 1) < size)
		{
		    os<<"mthread send to l + 1 "<<l + 1<<"\n";	
            destination = (root + l + 1) % size;
		    MPI_Isend(msb->send_buffer, msb->total_size, MPI_CHAR, destination, tag, MPI_COMM_WORLD, &msb->rrq);
			k++;
		}
        send_msg_count+=k;	
    } else {
        os<< "mthread send msg " << m->get_ray_size() <<"to "<<m->get_destination()<< std::endl;
        MPI_Isend(msb->send_buffer, msb->total_size, MPI_CHAR, m->get_destination(), tag, MPI_COMM_WORLD, &msb->lrq);
        k++;
	    send_ray_count++;
    }
    rs->accumulate_sent(k * m->get_ray_size());
	msb->n = k;
	mpi_in_flight.push_back(msb);
    
    return k;
}

void Communicator::collective(ProcStatus *rs) { }

bool Communicator::recv_message(RayList** List, RayStreamList * inList, ProcStatus *ps) {
    int recv_ready;
    MPI_Status status;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_ready, &status);

    if (recv_ready) {
        int count;
        MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &count);
        os<<"mthread recv msg "<<count<<"\n ";
        printf("%d before recv rays\n", rank);
        Message *recv_msg = new RecvMsg(List, inList, ps->get_local_chunk(), status, MPI_COMM_WORLD); 
        os<<"mthread recv rays from "<<recv_msg->get_sender()<<" size "<<recv_msg->get_ray_size() 
            <<" root: "<<recv_msg->get_root()<<"chunk "<<recv_msg->get_chunk()<<"\n";
        
        if (recv_msg->is_broadcast()) {
            os << "mthread recv broadcast\n";
            Export(recv_msg, ps); 
        }
       // os << "mthread switch msg type "<<recv_msg->get_type()<<"\n";
        if(isMaster()) 
            printf("master comm before recv\n");
        switch(recv_msg->get_type()) {
            case Quit: 
                { 
                    ps->set_exit();
                    return true;
                }
            case Status: 
                {
                    if(isMaster()) 
                        printf("master recv status\n");
                    int sender = std::max(recv_msg->get_sender(), recv_msg->get_root()); 
                        os << "mthread| recv status set proc " <<sender << "idle \n";
                    int * tmp = ps->get_status();
                        os <<tmp[0]<<" "<<tmp[1] <<" "<<tmp[2]<<" "<<tmp[3]<<"\n";
                    if(List[recv_msg->get_chunk()]->empty())
                        ps->set_proc_idle(sender);
                    
                    bool res = ps->update_global_rays((int*)recv_msg->get_content());
                    int *s = ps->get_status();
                    printf("recv status %d %d %d %d\n", s[0], s[1], s[2], s[3]);
                    os<< "update all rays received "<< s[0] <<" "<< s[1] <<" "<<s[2]<<" "<<s[3]<<"\n";
                    return true;
                } 
            case Ray: 
                {
                    if(isMaster()) 
                        printf("master recv ray\n");
                    os << "mthread| ray recv \n";
                    
                    ps->accumulate_recv(recv_msg->get_ray_size());
                    //set itself busy
                    ps->set_proc_busy(rank);
                    recv_ray_count++;
                    return true; 
                }
            case Schedule:
                {
                //    ps->set_exit();
                   // //if this proc need new chunk
                    os<<"recv schedule\n";
                    int * a = ps->get_chunk_map();
                    os<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<"\n";
                    if(ps->update_chunk((int*)recv_msg->get_content())) {
                        for(int i = 0; i < ps->get_chunk_size(); i++) {
                            List[i]->type = "out"; 
                        }
                        List[ps->get_local_chunk()]->type = "in";
                        os<<"mthread new chun k" << ps->get_local_chunk()<<"\n";
                        int * s = (int*)recv_msg->get_content();
                        os<<s[0]<<" "<<s[1]<<" "<<s[2]<<" "<<s[3]<<"\n";
                        int * a = ps->get_chunk_map();
                        os<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<"\n";
                        printf("mthread new chunk %d\n", ps->get_local_chunk());
                        
                    //    Message *recv_ray_msg = new RecvMsg(List, status, MPI_COMM_WORLD); 
                    //    ps->accumulate_recv(recv_msg->get_ray_size());
                        
                        return false; // waiting for rays
                    } else {
                        os<<"do nothing\n";
                        int * a = ps->get_chunk_map();
                        os<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<"\n";
                        return true;
                    }
                }
            default : 
                {
                    if(isMaster()) 
                        printf("master recv unknown\n");
                    os<<"recv unknown msg\n";
                    return true;
                }
        }
    } else {
        return false;
    } 
}

void Communicator::send_message(Message *outgoing_message, ProcStatus *rs) {

    //outList lock
    //serialize
    
    //current process is not destination or broadcast, send to destination 
    os << "mthread| send message"<<outgoing_message->get_destination() <<"\n";
    if (outgoing_message->is_broadcast() || (outgoing_message->get_destination() != rank))
        Export(outgoing_message, rs);

    os << "mthread| after send message\n";    
    if (outgoing_message->is_broadcast()) {
        if (outgoing_message->is_collective()) {
            rs->set_exit();
            os << "set exit\n";    
        } else {
            os << "mthread| broadcast schedule or status";    
        } 
    }
}

