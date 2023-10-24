//!
//! queue.h
//!
//! CS 121.Bolden
//! Clang 16.0.6
//! Caleb Barger
//! 10/20/23 
//! Linux x86_64
//! barg8397@vandals.uidaho.edu
//!
//! ScUUFd Queue 
//!

#ifndef QUEUE_H

#include "cabarger_cs121_include.h"

#define QueueNodePrototype(T) \
  struct QueueNode##T{ \
    T val; \
    QueueNode##T* next; \
  }; 

#define QueuePrototype(T) \
  QueueNodePrototype(T) \
	struct Queue##T { \
    QueueNode##T* head; \
    QueueNode##T* free_list; \
	}; \
  Queue##T queue##T##Init(); \
  void queue##T##Enqueue(Queue##T* queue, Arena* arena, T val); \
  T queue##T##Dequeue(Queue##T* queue);  

#define QueueImpl(T) \
  Queue##T queue##T##Init() { \
    return Queue##T{.head = 0, .free_list = 0}; \
  } \
  \
  void queue##T##Enqueue(Queue##T* queue, Arena* arena, T val) { \
    QueueNode##T* new_node; \
    if (queue->free_list == 0) { \
      new_node = arenaAlloc(arena, QueueNode##T, 1); \
    } else { \
      new_node = queue->free_list; \
      queue->free_list = queue->free_list->next; \
    } \
    new_node->val = val; \
    new_node->next = queue->head; \
    queue->head = new_node; \
  } \
  \
  T queue##T##Dequeue(Queue##T* queue) { \
    QueueNode##T* curr = queue->head; \
    if (curr == 0) \
      InvalidPath; \
    QueueNode##T* prev = curr;  \
    while(curr->next != 0) { \
      prev = curr; \
      curr = curr->next; \
    } \
    if (prev == curr) \
      queue->head = 0; \
    \
    prev->next = 0; \
    curr->next = queue->free_list; \
    queue->free_list = curr; \
    \
    return curr->val; \
  }

#define QUEUE_H
#endif
