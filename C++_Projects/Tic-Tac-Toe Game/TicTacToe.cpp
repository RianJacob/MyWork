
// Problem number 18, chapter 7
#include <iostream>
using namespace std;

void valid_input(char a[][3], int, int);
void display(char a[][3], int);
bool check(char a[][3], int);

int main()
{

  char a[3][3] = {{'*', '*', '*'},
		 {'*', '*', '*'},
		 {'*', '*', '*'}};
  int row=3, P1=1, P2=2;
 
  for (int i=0; i<5; i++)
    {
     
	  cout<<"Here is the current chart status: \n";
	  display(a, row);

	  
	  cout<<"\nPlayer #"<<P1<<" turn: \n";
	  valid_input(a, row, P1);
          cout<<"\nHere is the current chart status: \n";
	  display(a,row);
	  if (check(a,row))
	  {
	    cout<<"Player #"<<P1<<" Won!\n";
	   break;
	  }
	  if (i==4)
	    {
	      cout<<"Tie between the two players\n";
	      break;
	    }
	    cout<<"\nPlayer #"<<P2<<" turn: \n";
	  valid_input(a, row, P2);
	  if (check(a, row))
	   {
	     //cout<<"Here is the current chart status: \n";
	     display(a,row);
	     cout<<"Player #"<<P2<<" Won!\n";
	     break;
	   }	
    }
  return 0;
}

  
void valid_input(char a[][3], int row, int P)
  {
    int r, c;
    
	  cout<<"Enter the row and column number: ";
	  cout<<"\nrow number: ";
	  cin>>r;
	  while (r<0 || r>2)
	    {
	      cout<<"Invalid row number\n";
	      cout<<"Please enter row number between 0 and 2: ";
              cin>>r;
	    }  
	  cout<<"column number: ";
	  cin>>c;

	  while(c<0 || c>2)
	    {
	      cout<<"Invalid column number\n";
	      cout<<"Please enter column number between 0 and 2: ";
	      cin >>c;
	    }

	  while (a[r][c]=='X' || a[r][c]=='O')
	    {
	      cout<<"This position is already taken!!!\n";
	      cout<<"Here is the current status of the chart: \n";
	      display(a, 3);
	      cout<<"\nEnter row and column index that is not already taken\n";
	      cout<<"row: ";
	      cin>>r;

	      while(r<0 || r>2)
		{
		  cout<<"Invalid row number\n";
		  cout<<"Enter row number: ";
		  cin>>r;
		}

	      cout<<"column number: ";
	      cin>>c;

	      while (c<0 || c>2)
		{
		  cout<<"Invalid column number\n";
		  cout<<"Enter column number: ";
		  cin>>c;
		}
	    }
	  if (a[r][c]=='*' && P==1)
	    a[r][c] = 'X';
	  else if (a[r][c]=='*' && P==2)
	    a[r][c] = 'O';
      
  }


void display(char a[][3], int row)
{
  for (int i=0; i<row; i++)
    {
      for (int j=0; j<3; j++)
	{
	cout<<a[i][j]<<"  ";
	}
      cout<<endl;
    }
}

bool check (char a[][3], int row)
{
  // first row;
  if ((a[0][0]=='X' && a[0][1]=='X' && a[0][2]=='X') || (a[0][0]=='O' && a[0][1]=='O' && a[0][2]=='O'))
    return true;
  // second row
  if ((a[1][0]=='X' && a[1][1]=='X' && a[1][2]=='X') || (a[1][0]=='O' && a[1][1]=='O' && a[1][2]=='O'))
    return true;
  // thirs row
  if ((a[2][0]=='X' && a[2][1]=='X' && a[2][2]=='X') || (a[2][0]=='O' && a[2][1]=='O' && a[2][2]=='O'))
    return true;

  //first column
if ((a[0][0]=='X' && a[1][0]=='X' && a[2][0]=='X') || (a[0][0]=='O' && a[1][0]=='O' && a[2][0]=='O'))
    return true;
  //second column
if ((a[0][1]=='X' && a[1][1]=='X' && a[2][1]=='X') || (a[0][1]=='O' && a[1][1]=='O' && a[2][1]=='O'))
    return true;
  // third column
if ((a[0][2]=='X' && a[1][2]=='X' && a[2][2]=='X') || (a[0][2]=='O' && a[1][2]=='O' && a[2][2]=='O'))
    return true;

  //diagonal
if ((a[0][0]=='X' && a[1][1]=='X' && a[2][2]=='X') || (a[0][0]=='O' && a[1][1]=='O' && a[2][2]=='O'))
    return true;
  //anti-diagonal
if ((a[0][2]=='X' && a[1][1]=='X' && a[2][0]=='X') || (a[0][2]=='O' && a[1][1]=='O' && a[2][0]=='O'))
    return true;

  else
    return false;
}
