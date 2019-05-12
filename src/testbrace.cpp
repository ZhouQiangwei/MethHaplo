#include <stdio.h>
#include <string.h>

void Testbrace(char ch, int & count, int numberline){
	if((ch == '}') && (count == 0))  //
	{
		printf("unmatched at line: %d\n", numberline);
		return;
	}
	else if((ch == '}') && (count !=0))
	{
		count--;
	}
	else if(ch == '{')
	{
		count++;
	}
	if(count == 0)
		return;//printf("matched\n");
	else
		return ; //printf("unmatched at line: %d\n", numberline);
}

int main(int argc, char* argv[])
{
    int count=0;
    char* filename = NULL;
    filename = argv[1];
    FILE* cppFILE = fopen(filename, "r");
    int BATBUF=5000;char s2t[BATBUF];
	int numberline=0;
    while(fgets(s2t,BATBUF,cppFILE)!=0){
		int i;numberline++;
		for(i=0; i< strlen(s2t); i++){
			if(s2t[i] == '{' || s2t[i] == '}'){
				Testbrace(s2t[i], count, numberline);
			}
		}
    }
    if(count == 1) printf("unmatched at line: %d\n", numberline);
    fclose(cppFILE);
    return 0;
}
