#include "mpi.h"
#include "Message.h"
#include "ProcStatus.h"
#include <fstream>

struct mpi_send_buffer {
// If p2p, use left
    MPI_Request lrq;
    MPI_Request rrq;

    int n;
    char *send_buffer;
    int total_size;
};

class Communicator {

public:
    std::ofstream os;

    Communicator();
    ~Communicator(); 
   
    int get_size() { return size; } 
    int get_rank() { return rank; }
    int isMaster(){return rank == master;} 
    int get_master() {return master;}

    void reduce(float* film, float *reduce_buffer, int pixel_num);

    void all_gather_float(float *a, float *res, int size); 

    // broadcast or p2p. return sent number 1, in p2p, 0, 1 or 2 in bcast
    int  Export(Message * m, ProcStatus *rs); 

    bool recv_message(ProcStatus *ps, bool block); 

    bool process_message(Message *recv_msg, ProcStatus *ps); 
    
    void send_message(Message* msg, ProcStatus *rs); 

    void purge_completed_mpi_buffers(); 

    bool out_list_empty(int);

    RayStreamList *inList;  
    RayList *outArrayList;
    RayStreamList *outStreamList;
private:
    MPI_Status  sta[3];
    MPI_Request req[3];
    MPI_Request lrq, rrq;

    MPI_Comm Client_Comm;
    int rank, size, master;
   
    int send_ray_count, recv_ray_count;
    int send_msg_count, recv_msg_count;
        
    std::vector<mpi_send_buffer*> mpi_in_flight;
};


Communicator::Communicator() {
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    master = size - 1;
    
    printf("comm rank %d \n", rank);
    os     = std::ofstream("out/proc_" + std::to_string(rank));
    send_ray_count = 0; recv_ray_count = 0;
    send_msg_count = 0; recv_msg_count = 0;
}
    
Communicator::~Communicator() {
    printf("delete communicator\n");
    std::cout<<"msg send "<< send_msg_count<<"recv "<<recv_msg_count
      <<"ray send "<< send_ray_count<<"recv "<<recv_ray_count<<"\n";
    MPI_Finalize();
}

void Communicator::all_gather_float(float *a, float *res, int size) {
    MPI_Allgather(a, 1, MPI_FLOAT, res, 1, MPI_FLOAT, MPI_COMM_WORLD);
}

void Communicator::reduce(float* film, float *reduce_buffer, int pixel_num){
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(film, reduce_buffer, pixel_num, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
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
    
//    os << "mthread root " << root<< " rank " << rank <<" 2*d + 1 "<<2*d+1<< std::endl;
	// Only export if its either P2P or broadcast and this isn't a leaf
	if (m->is_broadcast() && (2*d + 1) >= size) {
//        os << "mthread root | " << root<< " rank " << rank <<" 2*d + 1 "<<2*d+1<< std::endl;
        return 0;
    }

	int tag = (MPI_TAG_UB) ? t % MPI_TAG_UB : t % 65535;
	t++;

	// If its a broadcast message, choose up to two destinations based
	// on the broadcast root, the rank and the size.  Otherwise, just ship it.
    
	struct mpi_send_buffer *msb = new mpi_send_buffer;
   
    statistics.start("run => message_thread => send_message => export");
    msb->total_size = m->get_header_size() + (m->has_content() ? m->get_content_size() : 0);
    msb->send_buffer = (char *)malloc(msb->total_size);
    memcpy(msb->send_buffer, m->get_header(), m->get_header_size());
//    os<< "mthread msb->total size " << msb->total_size << " " <<m->get_size()<<"broadcast"<<m->is_broadcast()<<std::endl;
    if (m->has_content())
        memcpy(msb->send_buffer + m->get_header_size(), m->get_content(), m->get_content_size());
    statistics.end("run => message_thread => send_message => export");
     
    if (m->is_broadcast()) {
        os<<"mthread broadcast\n";
        int l = (2 * d) + 1;
		int destination = (root + l) % size;
		MPI_Isend(msb->send_buffer, msb->total_size, MPI_CHAR, destination, tag, MPI_COMM_WORLD, &msb->lrq);
        os <<"m>"<<destination<<" "<< m->get_tag() <<"\n";
		k++;

		if ((l + 1) < size)
		{
		    os<<"mthread export to l + 1 "<<l + 1<<" tag "<<tag<<"\n";	
            destination = (root + l + 1) % size;
		    MPI_Isend(msb->send_buffer, msb->total_size, MPI_CHAR, destination, tag, MPI_COMM_WORLD, &msb->rrq);
            os <<"m>"<<destination<<" "<< m->get_tag() <<"\n";
			k++;
		}
        send_msg_count+=k;	
    } else {
        os<< "mthread send msg " << m->get_ray_size() <<"to "<<m->get_destination()<<" tag "<<tag<< std::endl;
        MPI_Isend(msb->send_buffer, msb->total_size, MPI_CHAR, m->get_destination(), tag, MPI_COMM_WORLD, &msb->lrq);
        os <<"m>"<<m->get_destination()<<" "<< m->get_tag() <<"\n";
        k++;
	    send_ray_count++;
    }
    rs->accumulate_sent(k * m->get_ray_size());
	msb->n = k;
	mpi_in_flight.push_back(msb);
    
    return k;
}

bool Communicator::out_list_empty(int n) {
    if(outArrayList == NULL)
        return outStreamList[n].empty();
    else
        return outArrayList[n].empty();
}

bool Communicator::process_message(Message *recv_msg, ProcStatus *ps) {

//        os<<"m<"<<status.MPI_SOURCE/*recv_msg->get_sender()*/<<" "<<recv_msg->get_tag() <<"\n";
//        os<<"mthread recv rays from "<<recv_msg->get_sender()<<" tag "<<recv_msg->get_tag()<<" size "<<recv_msg->get_ray_size() 
//            <<" root: "<<recv_msg->get_root()<<"chunk "<<recv_msg->get_chunk()<<"\n";
    
    statistics.start("run => message_thread => recv_message => broadcast");
    if (recv_msg->is_broadcast()) {
//            os << "mthread recv broadcast\n";
        Export(recv_msg, ps); 
    }
    statistics.end("run => message_thread => recv_message => broadcast");
   // os << "mthread switch msg type "<<recv_msg->get_type()<<"\n";
    if(isMaster()) 
        printf("master comm before recv\n");
    switch(recv_msg->get_type()) {
        case Quit: { 
            ps->set_exit();
            return true;
        }
        case Status: {
            statistics.start("run => message_thread => recv_message => StatusMsg");
            if(isMaster()) 
                printf("master recv status\n");
            int sender = std::max(recv_msg->get_sender(), recv_msg->get_root()); 
            if(out_list_empty(recv_msg->get_chunk()))
                ps->set_proc_idle(sender);
            
            bool res = ps->update_global_rays((int*)recv_msg->get_content());
            int *s = ps->get_status();
            printf("recv status %d %d %d %d\n", s[0], s[1], s[2], s[3]);
//                    os<< "update all rays received "<< s[0] <<" "<< s[1] <<" "<<s[2]<<" "<<s[3]<<"\n";
            statistics.end("run => message_thread => recv_message => StatusMsg");
            return true;
        } 
        case ArrayRay: {
            if(isMaster()) 
                printf("master recv ray\n");
            
            ps->accumulate_recv(recv_msg->get_ray_size());
            //set itself busy
            ps->set_proc_busy(rank);
            recv_ray_count++;
            return true; 
        }
        case StreamRay: {
            ps->accumulate_recv(recv_msg->get_ray_size());
            //set itself busy
            ps->set_proc_busy(rank);
            recv_ray_count++;
            return true;
        }
        case Schedule: {
            if(ps->update_chunk((int*)recv_msg->get_content())) {
//                        os<<"mthread new chun k" << ps->get_current_chunk()<<"\n";
                int * s = (int*)recv_msg->get_content();
//                        os<<s[0]<<" "<<s[1]<<" "<<s[2]<<" "<<s[3]<<"\n";
                int * a = ps->get_chunk_proc();
//                        os<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<"\n";
                printf("mthread new chunk %d\n", ps->get_current_chunk());

                return false; // waiting for rays
            } else {
                os<<"do nothing\n";
                int * a = ps->get_chunk_proc();
                os<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<"\n";
                return true;
            }
        }
        default : {
            if(isMaster()) 
                printf("master recv unknown\n");
            os<<"recv unknown msg\n";
            return true;
        }
    }
}

bool Communicator::recv_message(ProcStatus *ps, bool block) {
    int recv_ready;
    MPI_Status status;
    statistics.start("run => message_thread => recv_message => probe");
    if(block) 
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    else
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_ready, &status);

    statistics.end("run => message_thread => recv_message => probe");

    if (recv_ready || block) {
        int count;
        MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &count);
        printf("%d before recv rays\n", rank);
        
        statistics.start("run => message_thread => recv_message => RecvMsg");
        
        Message *recv_msg = new RecvMsg(outArrayList, outStreamList, inList, ps->get_current_chunk(), status, MPI_COMM_WORLD); 

        statistics.end("run => message_thread => recv_message => RecvMsg");

        //statistics.end("run => message_thread => recv_message => handle msg");
        
        return process_message(recv_msg, ps);
        
        //statistics.end("run => message_thread => recv_message => handle msg");
    } else {
        return false;
    } 
}

void Communicator::send_message(Message *message, ProcStatus *rs) {

    //outList lock
    //serialize
    
    //current process is not destination or broadcast, send to destination 
    os <<"mthread| send message to "<<message->get_destination() <<" tag "<< message->get_tag() <<"\n";
//    os<< "mthread send message: size "<<message->get_size()<<" broadcast "<<message->is_broadcast()<<std::endl;
    if (message->is_broadcast() || (message->get_destination() != rank))
        Export(message, rs);

//    os << "mthread| after send message\n";    
    if (message->is_broadcast()) {
        if (message->is_collective()) {
            rs->set_exit();
            os << "set exit\n";    
        } else {
            os << "mthread| broadcast schedule or status";    
        } 
    }
}

