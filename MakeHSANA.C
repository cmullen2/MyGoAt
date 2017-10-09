{
  TString PWD = gSystem->Getenv("$PWD");
  gSystem->Exec("cd $HSANA");
  gROOT->LoadMacro("THSParticle.C++");
//  gROOT->LoadMacro("THSDataReader.C++");
   // gSystem->Exe("ln -s THSParticle_C.so libTHSParticle.so");
  gSystem->Exec(TString("cd ")+PWD);
  //Done forget to include $HSANA in LD_LIBRARY_PATH : setenv LD_LIBRARY_PATH "$LD_LIBRARY_PATH":"$HSANA"
  //WILL LIKELY NEED TO MAKE SYMBOLIC LINK MANUALLY!!
}
