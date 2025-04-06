#include "NumberList.h"
#include <iostream>
using namespace std;



// appendNode means adding node at the end of the list
void NumberList::appendNode(double num)
{
  NodeList *newNode;
  NodeList *nodePtr;

  newNode = new NodeList;   // create new node
  
  newNode->value = num;
  newNode->next = nullptr;  //new node will be the last node

  // if there is no node(s) in the list
  if (!head)
    head = newNode;

  else
    {
      nodePtr = head;

  while(nodePtr->next)
    {nodePtr = nodePtr->next;}

    nodePtr->next = newNode;
    }
}

void NumberList::displayList() const
{
  NodeList *nodePtr;
  nodePtr = head;

  while(nodePtr)
    {
      cout<<nodePtr->value<<endl;
      nodePtr = nodePtr->next;
    }
}

void NumberList:: insertNode(double num)
{
  NodeList *newNode;       // to create new node
  NodeList *nodePtr;       // to move through the list
  NodeList *previousNode;  // to point to the previous node

  newNode = new NodeList;
  newNode->value = num;
  
  // if there is no nodes in the list, then the new node will be the only node in the list
  if (!head)
    {
      head = newNode;
      newNode->next = nullptr;
    }

  else
    {
      nodePtr = head;
      previousNode = nullptr;   //initilization

      while(nodePtr!= nullptr && nodePtr->value < num)
	{
	  previousNode = nodePtr;
	  nodePtr = nodePtr->next;
	}

      if (previousNode == nullptr)
	{
	  head = newNode;
	  newNode->next = nodePtr;
	}

      else
	{
	  previousNode->next = newNode;
	  newNode->next = nodePtr;
	}
    }
}

void NumberList::deleteNode(double num)
{

  NodeList *nodePtr;
  NodeList *previousNode = nullptr;   //initialization

  // if there is no nodes in the list
  if (!head)
    return;

  nodePtr = head;
  
  if (head->value == num)
    {
      
      head = nodePtr->next;
      delete nodePtr;
    }

  else
    {
      while (nodePtr != nullptr && nodePtr->value!=num)
	{
	  previousNode = nodePtr;
	  nodePtr = nodePtr->next;
	}

      if (nodePtr)
	{
	  previousNode->next = nodePtr->next;
	  delete nodePtr;
	}
    }
}


void NumberList::reverseList()
{
  NodeList *previousNode = nullptr;
  NodeList *nodePtr;
  NodeList *nextPtr;

  nodePtr = head;

  while(nodePtr)
    {
      
      nextPtr = nodePtr->next;
      nodePtr->next = previousNode;
      previousNode = nodePtr;
      nodePtr = nextPtr;
    }

  head = previousNode;
}

      

void NumberList:: displayNodes(NumberList &list1, NumberList &list2)
{

  // this function display one element from list1, and one element from list2, withput repetiton

  // this function basically display one element from beginning of list, and one element from the end of the list
  // assusing list2 is in reversal order of list1

  NodeList *ptr1 = list1.head;
  NodeList *ptr2 = list2.head;

  int counter = 0;
  while(ptr1)
    {
      counter++;
      ptr1 = ptr1->next;
    }

  // return ptr1 to beginning of list1
  ptr1 = list1.head;
   
  
  int index = 0;
  if (counter%2==0)
    {
      while (index<counter/2)
	{
	  cout<<"List1: "<<ptr1->value<<endl;
	  cout<<"List2: "<<ptr2->value<<endl;
	  ptr1 = ptr1->next;
	  ptr2 = ptr2->next;
	  index++;
	}
    }
  
  
  else
    {
      while (index<=counter/2)
	{
	  cout<<"List1: "<<ptr1->value<<endl;
	  if (index==counter/2)
	    break;
	  cout<<"List2: "<<ptr2->value<<endl;
	  ptr1 = ptr1->next;
	  ptr2 = ptr2->next;
	  index++;
	}
    }
}
  
   
          
NumberList::~NumberList()
{
  NodeList *nodePtr;     // to move through the list
  NodeList *nextNode;    // to poin to the next node
  // position nodePtr at the head of the list
  nodePtr = head;

  // while nodePtr is not at the end of the list...
  while (nodePtr)
    {

      // save a pointer to the next node
      nextNode = nodePtr->next;

      // delete the current node.
      delete nodePtr;

      // position nodePtr at the next node.
      nodePtr = nextNode;
    }
}


  
