#ifndef gnuplotH
#define gnuplotH

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>


using namespace std;

class gnuplot {

public:
    gnuplot(string gpfile="",bool pipe=false);
    ~gnuplot() {};
    
    void setterm(string terminfo);
    void addcommand(string cmd);
    
    void setxylable(string xlabel="x",string ylabel="y");
    void plotfile(string file,string parameters="",string range="");
    void replotfile(string file,string parameters="",string range="");
    
    void setout(string outfn="");
    
    void run(string cmdex="");
    void show()
    {addcommand("pause mouse close");run();};
    
    void reset() {cmds="";term="";};
    
private:
    bool usepipe;
    FILE *gpipe;
    string gpfn;
    string term;
    string cmds;
};

//-- end unit ----------------------------------------------------------------
#endif

