#include "NumberList.h"
#include <iostream>

using namespace std;

int main()
{
  NumberList list1;
  cout<<"Initially, List1 is empty.\n\n";

  // appending nodes, adding nodes to the end of the list
  cout<<"Now, appending new values to List1: 2.3, 5.6 "<<endl;
  list1.appendNode(2.3);
  list1.appendNode(5.6);
  cout<<"List1: \n";

  list1.displayList();
  cout<<endl;
  // the new nodes will be inserted in ascending order 
  cout<<"Inserting other nodes to List1: 1.1, 6.4, 5.5, 0.9 \n";
  cout<<"List1: "<<endl;
  list1.insertNode(1.1);
  list1.insertNode(6.4);
  list1.insertNode(5.5);
  list1.insertNode(0.9);
  list1.displayList();
  cout<<endl;
  
  
  //creating another identical list
     NumberList list2(list1);
     cout<<"List2 is identical copy of List1 "<<endl;
     cout<<"List2: "<<endl;
     list2.displayList();
  // reversing the order of elements in list2  
  list2.reverseList();
  cout<<endl;
  cout<<"List2 reversed: "<<endl;
  list2.displayList();
  cout<<endl;
  cout<<"List1 not changed "<<endl;
  cout<<"List1: "<<endl;
  list1.displayList();
  cout<<endl;
  
  cout<<"Now diplaying one element from list1, and one element from list2, without repetition"<<endl;
  list1.displayNodes(list1, list2); 
  
  cout<<endl;
  cout<<"now deleteing the node with value 5.5 from List1 \n";
  list1.deleteNode(5.5);
  cout<<"Here is the List1: \n";
  list1.displayList();
  cout<<endl; 

  cout<<"Now deleteing the node with value 0.9 from List1 "<<endl;
  list1.deleteNode(0.9);
  cout<<"Here is List1: "<<endl;
  list1.displayList();
  cout<<endl;
  cout<<"Now deleting the node with value 6.4 from List1 "<<endl;
  list1.deleteNode(6.4);
  cout<<"Here is List1: "<<endl;
  list1.displayList();
  
return 0;
}
