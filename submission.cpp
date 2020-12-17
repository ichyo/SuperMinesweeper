#include <iostream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

class SuperMinesweeper
{
  public:
    vector< vector<int> > Grid;
    int lastR;
    int lastC;

    SuperMinesweeper(int N, int M, int D, int row, int col)
    {
      Grid.resize(N, vector<int>(N));
      for (int r=0; r<N; r++)
        for (int c=0; c<N; c++)
          Grid[r][c]=-1;
          
      Grid[row][col]=0;

      lastR=0;
      lastC=-1;
    }

    string move()
    {        
      int r=lastR;
      int c=lastC;
      
      while (true)
      {
        c++;
        if (c>=Grid[r].size())
        {
          c=0;
          r++;
        }
        
        if (Grid[r][c]==-1) break;
      }
      
      lastR=r;
      lastC=c;
      return "G "+to_string(r)+" "+to_string(c);
    }

    void update(string feedback)
    {
      if (feedback.length()==0) return;
      
      int value;
      long runtime;
      stringstream ss(feedback);
      ss >> value >> runtime;
      
      Grid[lastR][lastC]=value;
    }
};

int main()
{
  string feedback;
  int N, M, D, row, col;
  cin >> N >> M >> D >> row >> col;
  getline(cin, feedback); // read endline

  SuperMinesweeper prog(N, M, D, row, col);
    
  while(true)
  {
    string move = prog.move();
    cout << move << endl;
    cout.flush();
    
    getline(cin, feedback);
    
    //stop if we have hit a mine
    if (feedback.find("BOOM!")!=string::npos)
    {
      cout << "STOP" << endl;
      cout.flush();
      break;
    }
    else
    {
      prog.update(feedback);
    }
  }
  return 0;
}
