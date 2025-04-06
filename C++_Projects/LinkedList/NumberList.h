#ifndef NumberList_H
#define NumberList_H
#include <iostream>
class NumberList
{
 private:
  struct NodeList
  {
    double value;
    NodeList *next;
  };

  NodeList *head;

 public:
  // Constructor 
  NumberList()
    { head = nullptr; }

  // Copy constructor 
  NumberList(NumberList &list)
  {
    head = nullptr;
    NodeList *ptr;
    ptr = list.head;
    
    while (ptr)
      {
	appendNode(ptr->value);
	ptr = ptr->next;
      }
  }
    
  // Destructor 
  ~NumberList();
    

  void appendNode(double);
  void displayList() const;
  void insertNode(double);
  void deleteNode(double);
  void reverseList();
  void displayNodes(NumberList &, NumberList &);  
};

#endif
