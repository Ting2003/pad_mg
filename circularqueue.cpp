#include "circularqueue.h"
using namespace std;

CircularQueue::CircularQueue(int maxSize):
	 queueHead_(0), queueTail_(0),
	 maxSize_(maxSize) {
         	queue_.resize(maxSize) ;
    	}

void CircularQueue::insert(Node* node) {
	if (queueTail_ >= queue_.size()) {
		cout << "Error in queue implementation" << endl ;
		exit(0) ;
	}
	queue_[queueTail_] = node ;
	queueTail_ = (queueTail_ + 1) % maxSize_ ;
	if (queueHead_ == queueTail_) {
		cout << "Queue is full. Increase maxSize of queue" 
		<< endl ;
		exit(0) ;
	}
}

Node* CircularQueue::extractFront(){
       if (queueHead_ == queueTail_) {
            cout << "Queue is empty. You are extracting more \
		elements than you inserted" <<  endl ;
           exit(0) ;
      }

       Node* x = queue_[queueHead_] ;
       queueHead_ = (queueHead_ + 1) % maxSize_ ;
       return x ;
}

bool CircularQueue::isEmpty(){
    return queueHead_ == queueTail_ ;
}

CircularQueue::~CircularQueue(){
     queueHead_ = 0 ;
     queueTail_ = 0 ;
     queue_.clear();
}

void CircularQueue::reset(){
     queueHead_ = 0 ;
     queueTail_ = 0 ;
}

