#include <cstring>
#include <cassert>
#include <netdb.h>
#include <iostream>
#include <sstream>
#include "trace_mq.h"

TraceMQ::TraceMQ(
  std::string r_dest_ip, int r_dest_port, 
  std::string r_controller_ip, int r_controller_port,
  int l_publisher_port,
  int comm_rank, int comm_size) : 
    r_dest_ip_ {r_dest_ip}, 
    r_dest_port_ {r_dest_port}, 
    r_controller_ip_ {r_controller_ip},
    r_controller_port_ {r_controller_port},
    l_publisher_port_ {l_publisher_port},
    comm_rank_ {comm_rank},
    comm_size_ {comm_size},
    state_ {TMQ_State::DATA},  /// Initial state is expecting DATA
    seq_ {0}
{
  std::string addrc("tcp://" + r_dest_ip_ + ":" + 
    std::to_string(static_cast<long long>(r_dest_port_+comm_rank_)));
  std::cout << "[" << comm_rank_ << "] Destination address=" << addrc << std::endl;

  context = zmq_ctx_new();
  /// setup the socket for data acquisition process
  server = zmq_socket(context, ZMQ_REQ); 
  zmq_connect(server, addrc.c_str()); 

  /// Setup sockets with controller process
  // PUB
  std::string addr("tcp://*:" + std::to_string(static_cast<long long>(l_publisher_port_)));
  std::cout << "[" << comm_rank_ << "] Local publisher address=" << addr << std::endl;
  publisher = zmq_socket (context, ZMQ_PUB);
  int rc = zmq_bind (publisher, addr.c_str()); assert (rc == 0);

  // REQ
  std::string addrcc("tcp://" + r_controller_ip_ + ":" + 
    std::to_string(static_cast<long long>(r_controller_port_)));
  std::cout << "[" << comm_rank_ << "] Destination controller address=" << addrcc << std::endl;
  controller = zmq_socket(context, ZMQ_REQ); 
  zmq_connect(controller, addrcc.c_str()); 
}

char* TraceMQ::MyHostname(){
  char *hostname = (char*) malloc(1024);
  memset((void*)hostname, '\0', 1024);
  gethostname(hostname, 1023);
  struct hostent* h;
  h = gethostbyname(hostname);
  memset((void*)hostname, '\0', 1024);
  /// h->h_name has the full hostname address
  sprintf(hostname, "%s", h->h_name);
  printf("My hostname is %s\n", hostname);

  return hostname;
}

void TraceMQ::Initialize() {
  /// Handshake with server
  tomo_msg_t *msg = prepare_data_info_req_msg(seq_, comm_rank_, comm_size_);
  send_msg(server, msg);
  free_msg(msg);
  ++seq_;

  /// Receive data info msg nad set local metadata structure
  msg = recv_msg(server); assert(seq_==msg->seq_n);
  metadata(*(read_data_info_rep(msg)));
  // release allocations
  free_msg(msg);
  ++seq_;

  /// Check if server has any projection
  msg = prepare_data_req_msg(seq_);
  send_msg(server, msg);
  free_msg(msg);
  ++seq_;

  /// Handshake with controller
  /// Send: hostname;rank;beg_sin_id;#sinogs;#cols
  /// use metadata above to setup message
  int m_size=2048;
  char *message= (char*) malloc(m_size);   /// FIXME:Constant message size
  memset((void*)message, '\0', m_size);
  char *hname = MyHostname();
  sprintf(message, "%s;%u;%u;%u;%u;%u;%u;%u;", hname, 
          l_publisher_port_, comm_rank_, comm_size_, 
          metadata().beg_sinogram, metadata().n_sinograms, 
          metadata().tn_sinograms, metadata().n_rays_per_proj_row);
  free(hname);
  zmq_send(controller, message, m_size, 0);
  free(message);
  printf("synching with controller..");
  /// TODO: At this point make sure controller connects to publisher of 
  /// this worker and then sends synch message.
  char tmp;
  int rc = zmq_recv(controller, &tmp, 0, 0); assert (rc == 0);  /// Synch with controller
  printf("done\n");
}


tomo_msg_t* TraceMQ::ReceiveMsg() {
  /// If previously fin message was recevied, return nullptr
  if(state()==TMQ_State::FIN) return nullptr;

  tomo_msg_t *dmsg = recv_msg(server);
  assert(seq_==dmsg->seq_n); ++seq_;
  if(dmsg->type == TRACEMQ_MSG_DATA_REP) { /// Message has data
    /// Tell data acquisition machine that you received the projection data
    tomo_msg_t *msg = prepare_data_req_msg(seq_);
    send_msg(server, msg);
    free_msg(msg);
    ++seq_;

    state(TMQ_State::DATA);
    return dmsg;
  } 
  else if (dmsg->type == TRACEMQ_MSG_FIN_REP){ /// Message is for finalization
    free_msg(dmsg);

    // send fin received msg
    tomo_msg_t *msg = (tomo_msg_t *) malloc(sizeof(tomo_msg_t));
    msg->type=TRACEMQ_MSG_FIN_REP;
    msg->seq_n=seq_;
    msg->size=sizeof(tomo_msg_t);
    send_msg(server, msg);
    free_msg(msg);
    ++seq_;

    state(TMQ_State::FIN);
    return nullptr;
  }
  std::cerr << "Unknown TRACEMQ_MSG state!" << std::endl;
  exit(1);
}

TraceMQ::~TraceMQ() {
  zmq_close(controller);
  zmq_close(publisher);
  zmq_close(server);
  zmq_ctx_destroy(context);
}

void TraceMQ::setup_msg_header(
  tomo_msg_t *msg_h, 
  uint64_t seq_n, uint64_t type, uint64_t size)
{
  msg_h->seq_n = seq_n;
  msg_h->type = type;
  msg_h->size = size;
}

tomo_msg_t* TraceMQ::prepare_data_req_msg(uint64_t seq_n)
{
  size_t tot_msg_size = sizeof(tomo_msg_t);
  tomo_msg_t *msg = (tomo_msg_t *) malloc(tot_msg_size);
  setup_msg_header(msg, seq_n, TRACEMQ_MSG_DATA_REQ, tot_msg_size);

  return msg;
}

tomo_msg_t* TraceMQ::prepare_data_rep_msg(uint64_t seq_n,
                                         int projection_id,        
                                         float theta,              
                                         float center,
                                         uint64_t data_size,
                                         float *data)
{
  uint64_t tot_msg_size=sizeof(tomo_msg_t)+sizeof(tomo_msg_data_t)+data_size;
  tomo_msg_t *msg_h = (tomo_msg_t *) malloc(tot_msg_size);
  tomo_msg_data_t *msg = (tomo_msg_data_t *) msg_h->data;
  setup_msg_header(msg_h, seq_n, TRACEMQ_MSG_DATA_REP, tot_msg_size);

  msg->projection_id = projection_id;
  msg->theta = theta;
  msg->center = center;
  memcpy(msg->data, data, data_size);

  return msg_h;
}

tomo_msg_data_t* TraceMQ::read_data(tomo_msg_t *msg){
  return (tomo_msg_data_t *) msg->data;
}

void TraceMQ::print_data(tomo_msg_data_t *msg, size_t data_count){
  printf("projection_id=%u; theta=%f; center=%f\n", 
    msg->projection_id, msg->theta, msg->center);
  for(size_t i=0; i<data_count; ++i)
    printf("%f ", msg->data[i]);
  printf("\n");
}



tomo_msg_t* TraceMQ::prepare_data_info_rep_msg(uint64_t seq_n, 
                                              int beg_sinogram, int n_sinograms,
                                              int n_rays_per_proj_row,
                                              uint64_t tn_sinograms)
{
  uint64_t tot_msg_size = sizeof(tomo_msg_t)+sizeof(tomo_msg_data_info_rep_t);
  tomo_msg_t *msg = (tomo_msg_t *) malloc(tot_msg_size);
  setup_msg_header(msg, seq_n, TRACEMQ_MSG_DATAINFO_REP, tot_msg_size);

  tomo_msg_data_info_rep_t *info = (tomo_msg_data_info_rep_t *) msg->data;
  info->tn_sinograms = tn_sinograms;
  info->beg_sinogram = beg_sinogram;
  info->n_sinograms = n_sinograms;
  info->n_rays_per_proj_row = n_rays_per_proj_row;

  return msg;
}
tomo_msg_data_info_rep_t* TraceMQ::read_data_info_rep(tomo_msg_t *msg){
  return (tomo_msg_data_info_rep_t *) msg->data;
}
void TraceMQ::print_data_info_rep_msg(tomo_msg_data_info_rep_t *msg){
  printf("Total # sinograms=%u; Beginning sinogram id=%u;"
          "# assigned sinograms=%u; # rays per projection row=%u\n", 
          msg->tn_sinograms, msg->beg_sinogram, msg->n_sinograms, 
          msg->n_rays_per_proj_row);
}

tomo_msg_t* TraceMQ::prepare_data_info_req_msg(uint64_t seq_n, 
                                               uint32_t comm_rank, 
                                               uint32_t comm_size)
{
  uint64_t tot_msg_size = sizeof(tomo_msg_t)+sizeof(tomo_msg_data_info_req_t);
  //printf("total message size=%llu\n", tot_msg_size);
  tomo_msg_t *msg = (tomo_msg_t *) malloc(tot_msg_size);
  setup_msg_header(msg, seq_n, TRACEMQ_MSG_DATAINFO_REQ, tot_msg_size);

  tomo_msg_data_info_req_t *info = (tomo_msg_data_info_req_t *) msg->data;
  info->comm_rank = comm_rank;
  info->comm_size = comm_size;

  return msg;
}
tomo_msg_data_info_req_t* TraceMQ::read_data_info_req(tomo_msg_t *msg){
  return (tomo_msg_data_info_req_t *) msg->data;
}

tomo_msg_t* TraceMQ::prepare_fin_msg(uint64_t seq_n){
  uint64_t tot_msg_size = sizeof(tomo_msg_t);
  tomo_msg_t *msg = (tomo_msg_t *) malloc(tot_msg_size);
  setup_msg_header(msg, seq_n, TRACEMQ_MSG_FIN_REP, tot_msg_size);

  return msg;
}

void TraceMQ::send_msg(void *server, tomo_msg_t* msg){
  zmq_msg_t zmsg;
  int rc = zmq_msg_init_size(&zmsg, msg->size); assert(rc==0);
  memcpy((void*)zmq_msg_data(&zmsg), (void*)msg, msg->size);
  rc = zmq_msg_send(&zmsg, server, 0); assert(rc==(int)msg->size);
}

tomo_msg_t* TraceMQ::recv_msg(void *server){
  zmq_msg_t zmsg;
  int rc = zmq_msg_init(&zmsg); assert(rc==0);
  rc = zmq_msg_recv(&zmsg, server, 0); assert(rc!=-1);
  /// Message size and calculated total message size needst to be the same
  /// FIXME?: We put tomo_msg_t.size to calculate zmq message size before it is
  /// being sent. It is being only being used for sanity check at the receiver
  /// side. 
  //printf("zmq_msg_size(&zmsg)=%zu; ((tomo_msg_t*)&zmsg)->size=%zu", zmq_msg_size(&zmsg), ((tomo_msg_t*)&zmsg)->size);
  //assert(zmq_msg_size(&zmsg)==((tomo_msg_t*)&zmsg)->size);

  tomo_msg_t *msg = (tomo_msg_t *) malloc(((tomo_msg_t*)zmq_msg_data(&zmsg))->size);
  /// Zero-copy would have been better
  memcpy(msg, zmq_msg_data(&zmsg), zmq_msg_size(&zmsg));
  zmq_msg_close(&zmsg);

  return msg;
}

void TraceMQ::free_msg(tomo_msg_t *msg) {
  free(msg); msg=NULL;
}


void TraceMQ::publish(float *data, int count){
  /// Prepare envelop 
  //zmq_msg_t zmsg_env;
  std::ostringstream label;
  label << metadata().beg_sinogram;
  //int rc = zmq_msg_init_size(&zmsg_env, label.str().size()); assert(rc==0);
  //memcpy( (void*)zmq_msg_data(&zmsg_env), 
  //        (void*)(label.str().data()), label.str().size());
  //rc = zmq_send(publisher, &zmsg_env, label.str().size(), ZMQ_SNDMORE);
  int rc = zmq_send(publisher, label.str().c_str(), label.length(), ZMQ_SNDMORE);
  std::cout << "Sending img size=" << count << std::endl;
  std::cout << "label size=" << label.str().length() << std::endl;

  /// Prepare and send the data with envelop
  zmq_msg_t zmsg;
  size_t data_size = sizeof(*data)*count;
  rc = zmq_msg_init_size(&zmsg, count*sizeof(float)); assert(rc==0);
  memcpy((void*)zmq_msg_data(&zmsg), (void*)data, data_size);
  //rc = zmq_msg_send(&zmsg, publisher, ZMQ_DONTWAIT); assert(rc==(int)data_size);  /// Non-blocking
  rc = zmq_msg_send(&zmsg, publisher, 0); assert(rc==(int)data_size); /// Blocking
  //rc = zmq_send(publisher, data, data_size, 0);
}
