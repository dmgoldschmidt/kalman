#include "GetOpt.h"

GetOpt::GetOpt(int argc, char** argv, char const* help_msg){
  bool opt_read = false; // parser state

  for(int i = 1;i < argc;i++){
    if(*argv[i] == '-' && *(argv[i]+1) != '\0'){
      argv[i]++; // skip the option flag
      if(strcmp("help",argv[i]) == 0){ // "-help" prints the help_msg and exits
        cout << help_msg;
        exit(0);
      }
      if(opt_read)argument.push_back(""); // two options in a row
      option.push_back(argv[i]);
      opt_read = true;
    }
    else{
      if(!opt_read)option.push_back(""); // two arguments in a row
      argument.push_back(argv[i]);
      opt_read = false;
    }
  }
  if(opt_read)argument.push_back(""); // we ended with an option
}

int GetOpt::get(char const* opt, int& arg){ // return integer argument
  int nmatches = get(opt);

  if(nmatches == 1) arg = atoi(argument[j].c_str());
  return nmatches;
}

int GetOpt::get(char const* opt, double& arg){ // return floating pt. argument
  int nmatches = get(opt);

  if(nmatches == 1) arg = atof(argument[j].c_str());
  return nmatches;
}

int GetOpt::get(char const* opt, float& arg){ // return floating pt. argument
  int nmatches = get(opt);

  if(nmatches == 1) arg = atof(argument[j].c_str());
  return nmatches;
}

int GetOpt::get(char const* opt, string& arg){// return string argument
  int nmatches = get(opt);

  if(nmatches == 1) arg = argument[j];
  return nmatches;
}

int GetOpt::get(char const* opt, char& arg){ // return char argument
  int nmatches = get(opt);
  
  if(nmatches == 1) arg = argument[j][0];
  return nmatches;
}

int GetOpt::get(char const* opt){ // look for option
  uint i,l;
  int nmatches = 0;

  for(i = 0;i < option.size();i++){
    l = strlen(option[i].c_str());
    if(l != 0 && strncmp(opt,option[i].c_str(),l) == 0){
      nmatches++;
      j = i;
    }
  }
  return nmatches;
}

  
  
#ifdef GetOpt_MAIN
int main(int argc, char**argv){
  GetOpt cmd_line(argc,argv, "testing help msg.\n");
  int a = cmd_line.argument.size();
  int o = cmd_line.option.size();
  int i_arg; // test integer argument (set for option = "int" below)
  double d_arg; // test double argument (option = "double")
  string s_arg; // test string argument (option = "string")
  // NOTE:  you can use any option names you want

  cout << o << " options:\n";
  for(int i = 0;i < o;i++)cout << cmd_line.option[i]<<endl;
  cout<<endl << a << " arguments:\n";
  for(int i = 0;i < a;i++)cout << cmd_line.argument[i]<<endl;
  if(cmd_line.get("int",i_arg) == 1) cout << "i_arg = "<< i_arg<<endl;
  if(cmd_line.get("double",d_arg) == 1) cout <<"d_arg = "<<d_arg<<endl;
  if(cmd_line.get("string",s_arg) == 1) cout <<"s_arg = "<<s_arg<<endl;
}
#endif
