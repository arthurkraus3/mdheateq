//---------------------------------------------------------------------------
// gnuplot.cpp
//---------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <string>

#pragma hdrstop

#include "fileutils.h"
#include "stringutils.h"

#include "gnuplot.h"

using namespace std;

//---------------------------------------------------------------------------
gnuplot::gnuplot(string gpfile,bool pipe) {
    usepipe=false; //not implemented
    if(gpfile=="") gpfn="gp.plt";
    else gpfn=gpfile;
    term="";
    cmds="";
};


void gnuplot::setterm(string terminfo) {
    term=terminfo;
};

void gnuplot::addcommand(string cmd) {
    cmds=cmds+cmd+"\n";
};

void gnuplot::setxylable(string xlabel,string ylabel) {
    if(xlabel!="") addcommand("set xlabel "+xlabel);
    if(ylabel!="") addcommand("set ylabel "+ylabel);
};

void gnuplot::plotfile(string file,string parameters,string range) {
    addcommand("plot "+range+" '"+file+"' "+parameters);
};

void gnuplot::replotfile(string file,string parameters,string range) {
    addcommand("replot "+range+" '"+file+"' "+parameters);
};

void gnuplot::setout(string outfn) {
    if(outfn!="") addcommand("set out '"+outfn+"'");
    else addcommand("set out");
};

void gnuplot::run(string cmdex) {
    fHandle f;
    string s;
    if(usepipe) return;
    
    s="";
    if(term!="") s="set term "+term+"\n";
    s=s+cmds;
    
    f=FileCreate(gpfn);
    FileWrite(f,s.c_str(),s.length());
    FileClose(f);
    
    s="gnuplot "+gpfn;
    if(cmdex!="") s=s+" "+cmdex;
    
    system(s.c_str());
};
