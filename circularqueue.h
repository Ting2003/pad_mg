#ifndef _CIRCULAR_QUEUE_H
#define _CIRCULAR_QUEUE_H
#include <vector>
#include <iostream>
#include <fstream>
#include "node.h"
class CircularQueue {
    int queueHead_, queueTail_ ;
    int maxSize_ ;
    vector<Node*> queue_ ;

public:
    CircularQueue(int maxSize); 
    ~CircularQueue();
    void insert(Node* node);       
    Node* extractFront(); 
    bool isEmpty();
    void reset(); 
} ;

#endif

