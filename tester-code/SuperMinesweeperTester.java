import java.awt.*;
import java.awt.geom.*;
import java.util.*;
import java.io.*;
import javax.imageio.*;

import com.topcoder.marathon.*;

public class SuperMinesweeperTester extends MarathonAnimatedVis {
    //parameter ranges
    private static final int minN = 10, maxN = 50;
    private static final double minM = 0.1, maxM = 0.3;
    private static final int minD = 1, maxD = 10;
    private static final int[] Dchoices = {1,2,4,5,8,9,10};

    //Inputs
    private int N;      //grid size
    private int M;      //number of mines    
    private int D;      //sensor distance
    private int zeroR;    //location of zero
    private int zeroC;  

    //Constants 
    private static final int MINE=-1;
    private static final int InitialReveals=1;

    //Graphics
    Image MinePic;
    Image FlagPic;

    //State Control
    private int[][] Grid;
    private boolean[][] Revealed;
    private boolean[][] Flagged;
    private boolean GameOver;    
    private int lastR;
    private int lastC;
    private int MaxNum;
    private int numFlagged;
    private int numGuesses;
    private int minesHit;



    protected void generate() {
        N = randomInt(minN, maxN);
        M = (int)(randomDouble(minM, maxM)*N*N);
        D = Dchoices[randomInt(0, Dchoices.length-1)];    

        //Special cases
        if (seed == 1)
        {
          N = minN;
          M = (int)(minM*N*N);
          D = 2;
        }
        else if (seed == 2)
        {
          N = maxN;
          M = (int)(maxM*N*N);
        }    
        

        //User defined parameters   
        if (parameters.isDefined("N"))
        {
          N = randomInt(parameters.getIntRange("N"), minN, maxN);
          M = (int)(randomDouble(minM, maxM)*N*N);
        }
        if (parameters.isDefined("M")) M = randomInt(parameters.getIntRange("M"), (int)(minM*N*N), (int)(maxM*N*N));        
        if (parameters.isDefined("D")) D = randomInt(parameters.getIntRange("D"), minD, maxD);

                
        java.util.List<Integer> zeros;        
        while(true)
        {
          Grid = new int[N][N];
          zeros = new ArrayList<Integer>();
          
          //place mines
          for (int placed=0; placed<M; )
          {
            int r=randomInt(0,N-1);
            int c=randomInt(0,N-1);
            if (Grid[r][c]!=MINE) 
            {
              Grid[r][c]=MINE;
              placed++;
            }
          }
          
          //compute number of mines around each empty cell
          for (int r=0; r<N; r++)
          {
            for (int c=0; c<N; c++)
            {
              if (Grid[r][c]==MINE) continue;
              
              for (int r2=0; r2<N; r2++)
                for (int c2=0; c2<N; c2++)
                  if (Grid[r2][c2]==MINE && sq(r2-r)+sq(c2-c)<=D)
                    Grid[r][c]++;
            }
          }
          
          //find zeros
          for (int r=0; r<N; r++)
            for (int c=0; c<N; c++)          
              if (Grid[r][c]==0)
                zeros.add(r*N+c);
                
          if (zeros.size()>0) break;
        }
        
        //display grid for debugging
        if (debug)
        {
          System.out.println("Grid:");          
          for (int r=0; r<N; r++)
          {
            for (int c=0; c<N; c++)        
              System.out.print((Grid[r][c]==MINE?"*":Grid[r][c])+" ");
            System.out.println();
          }        
        }
          
        
        //find max number
        MaxNum=0;
        for (int r=0; r<N; r++)
          for (int c=0; c<N; c++)
            if (Grid[r][c]!=MINE)
              MaxNum=Math.max(MaxNum,Grid[r][c]);
                
        
        Revealed = new boolean[N][N];        
        Flagged = new boolean[N][N];        
        GameOver = false;
        numFlagged = 0;
        numGuesses = InitialReveals;
        minesHit = 0;
        
        int ind = randomInt(0,zeros.size()-1);
        int zeroLoc = zeros.get(ind);
        zeroR = zeroLoc/N;
        zeroC = zeroLoc%N;
        Revealed[zeroR][zeroC]=true;
        
        lastR=-1;
        lastC=-1;
        

        if (debug) {
          System.out.println("Grid Size N = " + N);
          System.out.println("Mines M = " + M);
          System.out.println("Sensor Distance D = " + D);
          System.out.println("zero location: row "+zeroR+" column "+zeroC);
        }
    }
    
    private int sq(int a)
    {
      return a*a;
    }

    protected boolean isMaximize() {
        return true;
    }
    
    protected double run() throws Exception
    {
      init();
    
      if (parameters.isDefined("manual"))
      {
        setDefaultDelay(0);
        return 0;
      }
      else return runAuto();
    }
    
    
    protected double runAuto() throws Exception
    {
      double score = callSolution();
      if (score < 0) {
        if (!isReadActive()) return getErrorScore();
        return fatalError();
      }
      return score;    
    }    
    

    protected void contentClicked(double x, double y, int mouseButton, int clickCount)
    {
      if (!parameters.isDefined("manual")) return;
      if (GameOver) return;             //game is over
      
      int r=(int)Math.floor(y);
      int c=(int)Math.floor(x);
      if (!inGrid(r,c)) return;         //not in grid
            

      //flagging move - right click
      if (mouseButton == java.awt.event.MouseEvent.BUTTON3)
      {
        if (Flagged[r][c])
        {
          Revealed[r][c]=false;        
          Flagged[r][c]=false;
          numFlagged--;
        }
        else
        {
          if (Revealed[r][c]) return;
          
          Revealed[r][c]=true;
          Flagged[r][c]=true;
          numFlagged++;
        }
      }
      //guessing move - left click
      else if (mouseButton == java.awt.event.MouseEvent.BUTTON1)
      {        
        if (Revealed[r][c]) return;       //already revealed
        
        Revealed[r][c]=true;
        lastR=r;
        lastC=c;
        
        //we blew up
        if (Grid[r][c]==MINE)
        {
          minesHit++;
        }
        //we guessed a number
        else
        {
          numGuesses++;
        }
        
        double score=(numGuesses*1.0/(N*N - M))/(minesHit+1);
        addInfo("Score", String.format("%.3f", score)); 
      }
      
      addInfo("Flagged", numFlagged);
      addInfo("Guessed", numGuesses);
      addInfo("Mines Hit", minesHit);
      updateDelay();
    }    
    
    protected void timeout() {
        addInfo("Score", "-1 (Timeout)");
        addInfo("Time", getRunTime());
        update();
    }    
    
    protected boolean inGrid(int r, int c)
    {
      return (r>=0 && r<N && c>=0 && c<N);
    }
    
    private double callSolution() throws Exception {
        writeLine(N);
        writeLine(M);
        writeLine(D);
        writeLine(""+zeroR+" "+zeroC);
        flush();
        if (!isReadActive()) return -1;
        
        updateDelay();

        int numMoves = 0;
        double score = 0;        
        while ((numGuesses < N*N - M) && !GameOver)
        {
          startTime();
          String move = readLine();
          stopTime();
          
          if (move.equals("STOP"))
          {
            GameOver=true;
            if (debug)
            {
              System.out.println("Move #" + numMoves + ": " + move);
            }
            update(score);
            break;
          }
                 
          //get move
          int r,c;
          char type;
          try
          {
            String[] temp=move.split(" ");
            if (temp.length!=3)
            {
              setErrorMessage("Error reading move: "+move);
              return -1;            
            }
            type=temp[0].charAt(0);
            r=Integer.parseInt(temp[1]);
            c=Integer.parseInt(temp[2]);
          }
          catch(Exception e)
          {
            setErrorMessage("Error reading move: "+move);
            return -1;
          }
          if (!(type=='F' || type=='G'))
          {
            setErrorMessage("Unknown move type: "+move);
            return -1;          
          }
          if (!inGrid(r,c))
          {
            setErrorMessage("Invalid coordinates of move: "+move);
            return -1;          
          }
          if (Revealed[r][c])
          {
            setErrorMessage("Location of move has already been revealed: "+move);
            return -1;                    
          }
                                          
          
          //if we are here then the move was valid
          String feedback="";
          Revealed[r][c]=true;
          lastR=r;
          lastC=c;
          numMoves++;
          
          //flagging move
          if (type=='F')
          {
            Flagged[r][c]=true;
            numFlagged++;
          }
          //guessing move
          else
          {          
            //we blew up
            if (Grid[r][c]==MINE)
            {
              feedback="BOOM! "+getRunTime();
              minesHit++;
            }
            //we guessed a number
            else
            {
              feedback=Grid[r][c]+" "+getRunTime();
              numGuesses++;
            }
          }
          
          score=(numGuesses*1.0/(N*N - M))/(minesHit+1);
          update(score);

          if (debug)
          {
            System.out.println("Move #" + numMoves + ": " + move);
            System.out.println("Feedback #" + numMoves + ": " + feedback);
          }            

          writeLine(feedback);
          flush();
        }
        GameOver=true;

        return score;
    }
     
    private void update(double score) {
        if (hasVis()) {      
          synchronized (updateLock) {
            addInfo("Flagged", numFlagged);
            addInfo("Guessed", numGuesses);
            addInfo("Mines Hit", minesHit);
            addInfo("Score", String.format("%.3f", score));
            addInfo("Time", getRunTime());              
          }
          updateDelay();
        }
    }

    protected void paintContent(Graphics2D g) {
           
      adjustFont(g, Font.SANS_SERIF, Font.PLAIN, String.valueOf((""+MaxNum).length()), new Rectangle2D.Double(0, 0, 1, 1));    
      g.setStroke(new BasicStroke(0.005f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));              

      //draw grid
      for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++)
        { 
          g.setColor(Color.white);
          g.fillRect(c, r, 1, 1);
          
          if (Revealed[r][c])
          {
            if (GameOver && Flagged[r][c] && Grid[r][c]!=MINE)    //crossed out bomb
            {
              g.drawImage(MinePic,c,r,1,1,null);
              g.setColor(Color.red);
              GeneralPath gp = new GeneralPath();
              gp.moveTo(c+0.1,r+0.1); 
              gp.lineTo(c+0.9,r+0.9);
              gp.moveTo(c+0.9,r+0.1);
              gp.lineTo(c+0.1,r+0.9); 
              g.setStroke(new BasicStroke(0.15f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));   
              g.draw(gp);                
              g.setStroke(new BasicStroke(0.005f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));        
            }
            else if (Flagged[r][c])  g.drawImage(FlagPic,c,r,1,1,null);
            else if (Grid[r][c]==MINE)  g.drawImage(MinePic,c,r,1,1,null);
            else
            {
              g.setColor(Color.black);   
              drawString(g, String.valueOf(Grid[r][c]), new Rectangle2D.Double(c+0.5, r+0.5, 0, 0));            
            }
          }
                    
          g.setColor(Color.gray);                        
          g.drawRect(c, r, 1, 1);          
        }
        
      if (lastR!=-1)
      {
        g.setColor(Color.red);
        g.setStroke(new BasicStroke(0.05f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        g.drawRect(lastC, lastR, 1, 1);
      }
    }


    private void init() {
      if (hasVis())
      {
        MinePic = loadImage("images/bomb.png");
        FlagPic = loadImage("images/flag.png");
      
        setContentRect(0, 0, N, N);
        setInfoMaxDimension(15, 11);
        setDefaultDelay(200);

        addInfo("Seed", seed);
        addInfoBreak();
        addInfo("Size N", N);
        addInfo("Mines M", M);
        addInfo("Distance D", D);
        addInfoBreak();           
        addInfo("Flagged", "-");
        addInfo("Guessed", "-");
        addInfo("Mines Hit", "-");
        addInfo("Score", "-");
        addInfoBreak();
        addInfo("Time", "-");
        update();
      }
    }
    
    Image loadImage(String name) {
      try{
        Image im=ImageIO.read(new File(name));
        return im;
      } catch (Exception e) { 
        return null;  
      }             
    }       


    public static void main(String[] args) {
        new MarathonController().run(args);
    }
}
