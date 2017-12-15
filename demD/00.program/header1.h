/*
"fsimII.04.c"
FsimII
-------
Program created :	05/15/03	Anoop John
Last modified	:	07/21/03	Anoop John

Pseudocode
----------
{
	Read Control File
	Initialize data and counts
	Run Algorithm for initRunGen generations
		Repeat for repNum replicates
		{
			Run Algorithm for preRunGen generations
				Run Algorithm for repRunGen generations 
					Collect Data from replicate and output to file
					}
}
Algorithm
---------
{
	doMutations
	{
		For all individuals
		{
			For 1, 2, 3 codon positions
			{
				1->0 mutations	( ! remember to make same changes in both sections)
				0->1 mutations	( ! remember to make same changes in both sections)
			}
		}
	}
	doRecombination
	{
		For all pairs
		{
			find and sort positions for recombination
				swap alternate segments
			{
				if polymorphic site
					then if different
						then exchange data
						}
		}
	}
	getNextGenSeqs
	{
		Find Fitness of individuals
		Do multinomial sampling
		create next generation
	}
	updateCounts
	{
		Update counts of data being tracked
	}
}
*/
/*---------------------------------------------------------------------------------------*/
#define DEBUGMODE
#undef DEBUGMODE

#define DEBUGMODE_1					/* check mutation function and basects */
#undef DEBUGMODE_1
/**/
#define VERBOSE
#undef VERBOSE
/**/
/*---------------------------------------------------------------------------------------*/
/*constants*/
#define FLOATPRECISION	0.000001			/*precision to be used for calculations*/
#define PI				3.1415926535897932384626433832795
#define MAXLINECHAR		1000000				/*maximum characters in a line*/
#define	MAXFILENAME		501					/*maximum length of filename*/
#define IDUM_INIT		-123341				/*value used to initialize the random number generator*/

/*---------------------------------------------------------------------------------------*/
/*debug tools*/
#define pf(a) 		printf("%.20f\t", a)
#define pl(a) 		printf("%ld\t", a)
#define ps(a) 		printf("%s\t", a)
#define pa(a) 		printf("%c\t", a)
#define tab			printf("\t")
#define ent			printf("\n")
#define help		printf(" O'er h're")
#define halt		do{printf("\nPress enter to continue!\n");getchar();}while(0)
// #undef halt 
// #define halt		do{printf("\nPress any key to continue!\n");}while(0)

/*---------------------------------------------------------------------------------------*/
/*Error function implemented as macro*/
#define errorOut(a)	do{printf("\nError! ");printf a;printf("\nexiting program.....");getchar();exit(1);}while(0)
/*---------------------------------------------------------------------------------------*/
/*macro function to find codon position given site position*/
#define sitePosToCodPos(a)	((((a)+2)%3)+1)
/*---------------------------------------------------------------------------------------*/
/*integer random number generator*/
#define getRandLong_(min, max)	((min)+(long)((1.0-ran2())*((max)-(min)+1)))	/* careful about which ran2 is being used !!!!!!!!!!!!! */
/*---------------------------------------------------------------------------------------*/
long idumVal;						/*default ran2 seed. if seed given!=0 then it will be used*/	
long *idum = &idumVal;						

/*---------------------------------------------------------------------------------------*/
/*function declarations*/
double ran2(void);
double poidev(double xm);
double gammln(double xx);
void multdev(double inn[], long k, long n, long nn[]);
double bnldev(double pp, long n);

FILE* fileOpen(char *path, char *mode);
char* getFileString(FILE *fp, char *string, long maxStringLength);
long readDataWithComments(FILE* fp, char*format, void* target);
void* memAlloc(size_t nmemb, size_t membsize, char* membname);
char* getUserString(char *question, char *string, unsigned long maxStrLen);
long isBlankString(const char *string);
long getFileLineRemoveComments(FILE *fp, char *string, unsigned long maxLength);
char *extractStringFromQuotedString(char *string);
long strCmpIgnoreCase(const char *string1, const char *string2);
char *removeExtraWhiteSpaces(char *string);
char* getTokenNameAndValue(const char *string, char *token, const char *separator, char *value);
unsigned long programTimer(long m);
char * CreateFolder(const char *full_path, int mode_type);
void* memoryHandler(size_t memSize, void *orig, const char *variableName, char request);
long StringToLong(char * InputString);
double StringToDouble(char * InputString);
char* fixFolderName(char *folderName, const char *dirSep);
char * get_substr_str_tag_delim(char *parent_str_orig, char *tag, long *start_index, char *delimiter_str, long max_substring);
long count_char_in_str(char *input_string, char char_to_count);
long get_stringpos(char *parent, char *tag, long start_index);
char *get_line_from_file_rev( FILE *fptr, char *filename, long *eof_found);
long get_file_dataline_ct(char readfile[], long *last_line_eoln);

int compareUL(const void*a, const void*b);
long getTime(long);

/*---------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*safe memory allocation function. tries to allocate memory and if error then 
 print error with variable name and exit
 nmemb - the number of elements to be allocated
 memmsize - size of each element
 membname - name of the variable for which memory is being allocated for errorout
 */
void* memAlloc(size_t nmemb, size_t membsize, char* membname)
{
	void* ptr;
	static long totmem=0;
	if(!nmemb&&!membsize)
		return &totmem;
	ptr=calloc(nmemb, membsize);
	if(ptr)
	{
		totmem+=nmemb*membsize;
		return ptr;
	}
	else
		errorOut(("Unable to allocate %ld bytes for %s", nmemb*membsize, membname));
	return ptr;
}
/*---------------------------------------------------------------------------------------*/
/*function to free memory allocated through memAlloc or memRealloc. if pointers allocated
through calls to malloc, calloc or realloc are passed to this function then function will
print error and exit. It is advisable that all memory allocated through memAlloc and
memRealloc are freed through this function. otherwise the internal list in memoryHandler
will get corrupted.
ptr				- the pointer to be freed
created 		08-14-03	Anoop John
last Modified	08-14-03	Anoop John
from HALabUtils.h
*/
void freeMem(void *ptr)
{	
	memoryHandler(0, ptr, 0, 3);
	#ifdef VERBOSE
		memoryHandler(0, 0, 0, 4);
	#endif
	return;
}
/*--------------------------------------------------------------------------------*/
/*Safe file open function. can check for OS specific length issues with filenames
 for mac the length is 31 and for OSX it is 255
 path - the path to the file+filename or just filename if in same folder
 mode - the fopen modes
 */
FILE* fileOpen(char *path, char *mode)
{
	FILE *fp;
	long i, j;
	
	i=strlen(path);
	j=0;
	while(i>0)
	{
		i--;
		if(path[i]==':'||path[i]=="\\"[0]||path[i]=='/')
			break;
		j++;
	}
	if(j>MAXFILENAME)
		errorOut(("Filename %s exceeded the maximum limit set by the operating system", &path[i]));
	fp=fopen(path, mode);
	if(fp)
		return fp;
	else
		errorOut(("Unable to open %s in %s mode", path, mode));
	return fp;
}
/*---------------------------------------------------------------------------------------*/
/*Function that asks the user a question and gets a string including whitespaces
question - the question to be asked to the user
string - the char array where the string read is to be stored
maxStrLen - the maximum length of string
created 		07-20-03	Anoop John
last Modified	07-20-03	Anoop John
from HALabUtils.h
*/
char* getUserString(char *question, char *string, unsigned long maxStrLen)
{	
	char *origString;				/*pointer to store the original char pointer*/
	char c;							/*temporary character used to read char by char*/
	unsigned long strLength;		/*length of the string*/
	origString = string;
	if(string == NULL || question == NULL)
		errorOut(("Null Pointer passed to function"));
	printf("%s", question);
	strLength = 0;
	c=getchar();
	while(c != '\n' && c != '\r')
	{	
		strLength++;
		if(strLength >= maxStrLen)
			errorOut(("Exceeded the maximum length(%lu) of the string that can be read", maxStrLen));
		*string++ = c;
		c = getchar();
	}
	*string = 0;
	return origString;
}
/*--------------------------------------------------------------------------------*/
/*function to safely read a string with spaces from a file. checks for the maximum
 size of the string that can be read into the variable passed
 fp - file pointer at position where string is to be read
 string - the pointer to the location where the string is to be stored
 maxStringLength - the maximum length of the string that can be read into the variable
 created 		:06/06/03	- Anoop John 
 last modified	:06/06/03	- Anoop John
 */
char* getFileString(FILE *fp, char *string, long maxStringLength)
{
	long curLength;
	curLength=0;
	
	/*check for null pointer*/
	if(string==0)
		errorOut(("Null pointer passed to getFileString"));
	
	/*try to read till max length or end of file*/
	while(curLength<maxStringLength&&!feof(fp))
	{
		/*read char by character*/
		string[curLength]=fgetc(fp);
		/*read till end of line is read. break at end of line*/
		if(string[curLength]=='\n'||string[curLength]=='\r')
		{
			string[curLength]=0;
			break;
		}
		else
		{
			/*if new line is not reached after max length has been read then error*/
			curLength++;
			if(curLength==maxStringLength)
				errorOut(("String Length exceeded maximum allowed(%lu) in getFileString", maxStringLength));
		}
	}
	/*if end of file is reached error*/
	if(feof(fp))
		errorOut(("Unexpected end of file in getFileString"));
	return string;
}
/*---------------------------------------------------------------------------------------*/
/*checks if the null terminated string passed to the function has only whitespace characters
returns 1 if the string is empty or has only blank lines else returns 0
+ added error checks for strings [08/11/03 - Anoop]
string - the string to be checked
created 		07-01-03	Anoop John
last Modified	08-11-03	Anoop John
from HALabUtils.h
*/
long isBlankString(const char *string)
{
	/*added error check for null string [added 08/11/03 - Anoop]*/
	if(string == 0)
		errorOut(("Null Pointer passed to function"));
	/*isspace('\0') will return false*/
	/*read till null char or non whitespace char*/
	while(isspace(*string))
		string++;
	/*if end of string, string is blank return 1*/
	if(*string=='\0')
		return 1;
	/*if not end of string, string is not blank return 0 */
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*function to read a line from a file and remove c formatted comments at the end of the line
comments can span multiple lines. both types of C99 / * and / / commenting styles are permitted. 
string should be long enough to hold maxLength characters including the null char
+ added error checks for strings [08/11/03 - Anoop]
string - the char array into which the string read is to be stored
maxlength - the length of the string
Created 		06/06/03	Anoop John
Last Modified	08/11/03	Anoop John*/
long getFileLineRemoveComments(FILE *fp, char *string, unsigned long maxLength)
{	
	unsigned long curPos;
	char curChar;
	char lastChar;
	char commentStart;
	char commentEnd;

	/*added error check for null string [added 08/11/03 - Anoop]*/
	if(string == 0)
		errorOut(("Null Pointer passed to function"));
	/*to account for the null character*/
	maxLength--;
	commentStart = 0;
	commentEnd = 0;
	/*read a char from the file*/
	curChar = fgetc(fp);
	lastChar = 0;
	curPos = 0;
	/*read till commentEnd is found*/
	while(!feof(fp) && !commentEnd)
	{	
		/*check if comment has started commentStart == 1 means // type comment == 2 means * style comment*/
		if((lastChar == '/') && (curChar == '/'))
			commentStart = 1;
		if((lastChar == '/') && (curChar == '*'))
			commentStart = 2;
		/*If comment has not started then copy the char to the string*/
		if(commentStart == 0)
		{	
			/*check if length has exceeded maximum length*/
			if(curPos == maxLength)
				errorOut(("Line length exceeded the maximum in getFileLineRemoveComments"));
			/*if length is okay copy the current char into string*/
			string[curPos] = curChar;
			curPos++;
		}
		/*store current as the last char read*/
		lastChar = curChar;
		/*read next char*/
		curChar = fgetc(fp);
		/*if the comment had started check if it has ended*/
		/*if no comment and end of line or // comment and end of line then comment ended*/
		if(((commentStart == 1)||(commentStart == 0)) && ((curChar == '\n')||(curChar == '\r')))
		{
			commentEnd = 1;
		}
		/*if * type comment then check if found the ending */
		else if((commentStart == 2) && (curChar == '/') && (lastChar == '*'))
		{	
			commentEnd = 1;
			/*if ending is found then read whitespaces on the same line*/
			lastChar = curChar;
			curChar = fgetc(fp);
			
			while(curChar != '\n' && curChar != '\r' && !feof(fp))
			{
				lastChar = curChar;
				curChar = fgetc(fp);
				/*check to see if the character after the comment is non whitespace
				if so then error*/
				if(!isspace(curChar))
				{
					/*null terminate the string to print out properly*/
					if(commentStart == 0)
						string[curPos] = 0;
					else
						string[curPos-1] = 0;
					errorOut(("Non white-space character after comment after line %s", string));
				}
			}
		}		
	}
	/*null terminate the string*/
	if(commentStart == 0)
		string[curPos] = 0;
	else
		string[curPos-1] = 0;
	/*to support CRLF lineending if last char read is \r then try to see if next char is \n*/
	if(curChar == '\r')
	{
		lastChar = curChar;
		curChar = fgetc(fp);
		/*if not end of file and current character read is not a LF then go back in file*/
		if(!feof(fp) && curChar != '\n')
			fseek(fp, -1, SEEK_CUR);
	}
	/*if type 2 comment has started and ended before end of file then error*/
	if(commentStart == 2 && commentEnd == 0 && feof(fp) )
		errorOut(("Unexpected end of file before end of comment after line %s", string));
	return 0;
}
/*---------------------------------------------------------------------------------------*/
/*function  to  extract  a  string from a quoted string. the unquoted string will be placed
into  the  same variable. For example  if  string = "\"parameter\"", after  processing  the 
string = "parameter". So the quotes will be removed
string 			- the quoted string to be processed
created 		08-13-03	Anoop John
last Modified	08-13-03	Anoop John
from HALabUtil.h
*/
char *extractStringFromQuotedString(char *string)
{
	size_t curStrPos;					/*variable to move along the string*/
	size_t startPos;					/*the startPos of the stringValue*/
	size_t endPos;						/*the ending position of the stringValue*/
	size_t length;						/*length of the string*/
	
	if(string == 0)
		errorOut(("Null Pointer passed to function"));
	/*find the opening quote*/
	curStrPos = 0;
	while(string[curStrPos] != '\"')
	{
		if(isspace(string[curStrPos]) == 0 || string[curStrPos] == 0)
			errorOut(("%s is not a quoted string or it has non-white-space characters before opening quote", string));
		curStrPos++;
	}
	curStrPos++;
	startPos = curStrPos;
	/*find the ending quote*/
	while(string[curStrPos] != '\"')
	{	
		if(string[curStrPos] == 0)
			errorOut(("Unexpected end of string(%s) before closing quote", string));
		curStrPos++;
	}
	endPos = curStrPos;
	length = endPos-startPos;
	/*move characters*/
	memmove(string, &string[startPos], length*(sizeof(char)));
	/*null terminate the string*/
	string[length] = 0;
	return string;
}
/*---------------------------------------------------------------------------------------*/
/*function to compare two strings ignoring the case of alphabetical characters.
function returns a negative value if string1 is alphabetically on top and
positive value if string2 is alphabetically on top. if both strings are equal
the function returns 0. function assumes that the strings are null terminated.
+ added error checks for strings [08/11/03 - Anoop]
string1 - null terminated 0 based string
string2 - null terminated 0 based string
created 		06-27-03	Anoop John
last Modified	08-11-03	Anoop John
from "HALabUtils.h"
*/
long strCmpIgnoreCase(const char *string1, const char *string2)
{
	/*added error check for null string [added 08/11/03 - Anoop]*/
	if(string1 == 0 || string2 == 0)
		errorOut(("Null Pointer passed to function"));
	while((toupper(*string1) == toupper(*string2)) && *string1 && *string2)
	{
		string1++;
		string2++;
	}
	/*at this point both strings are different and one or both may be zero*/
	if(isalpha(*string1) && isalpha(*string2))
		return(toupper(*string1)-toupper(*string2));
	else
		return(*string1-*string2);
}
/*---------------------------------------------------------------------------------------*/
/*function to remove extra white spaces in a string. This function removes white 
spaces from the beginning and end of a string
returns the a pointer to the string as well
+ added error checks for strings [08/11/03 - Anoop]
+ added support for all whitespace strings [08/21/03 - Anoop]
string - the string from which the starting and ending whitespaces are to be removed
created 		07-01-03	Anoop John
last Modified	08-21-03	Anoop John
from "HALabUtils.h"
*/
char *removeExtraWhiteSpaces(char *string)
{
	char *newString;				/*the new string with newlines as LF*/
	unsigned long strStart;			/*the position of first nonwhitespace in the source string*/
	unsigned long strEnd;			/*the position of last nonwhitespace in the source string*/
	
	/*error check for null string [added 08/11/03 - Anoop]*/
	if(string == 0)
		errorOut(("Null Pointer passed to function"));
		
	#ifdef VERBOSE
		printf("before conversion :%s:\n", string);
	#endif
	newString = (char*) memAlloc (strlen(string)+1, sizeof(char), "newString :removeExtraWhiteSpaces");
	strStart = strEnd = 0;
	/*read until first nonwhitespace character - will read till end of string if string is all
	white space*/
	while(isspace(string[strStart]) != 0)
		strStart++;
	#ifdef VERBOSE
		printf("stringstart :%s:\n", &string[strStart]);
		printf("strstart %lu\n", strStart);
	#endif
	/*find end of string*/
	while(string[strEnd])
		strEnd++;
	/*move to the last char before null char*/
	/*added check for handling null strings [07/24/03*/
	if(strEnd > 0)
		strEnd--;
	#ifdef VERBOSE
		printf("strend %lu\n", strEnd);
	#endif
	/*go back and find first nonwhitespace character*/
	/*added check for handling null strings [07/24/03 Anoop]*/
	while(isspace(string[strEnd]) != 0 && strEnd > 0)
		strEnd--;
	#ifdef VERBOSE
		printf("strend %lu\n", strEnd);
		printf("stringEnd :%s:\n", &string[strEnd]);
	#endif
	/*if the string is full whitespace return an empty string "" [added 08/21/03 Anoop]*/
	if(string[strStart] == '\0' && strEnd == 0)
	{	
		newString[0] = '\0';
	}
	else
	{
		/*strStart and strEnd points to first and last nonwhitespace*/
		strncpy(newString, &string[strStart], strEnd-strStart+1);
		/*strncpy doesnt null terminate if all n chars are read from source*/
		newString[strEnd-strStart+1] = 0;
	}
	strcpy(string, newString);
	free(newString);
	#ifdef VERBOSE
		printf("after conversion :%s:\n", string);
	#endif
	return string;
}
/*--------------------------------------------------------------------------------*/
/*function to extract token names and values from a string with token on the left and
values on the right. function to be used to read control files. anything on the left of
the first separator will be taken as the token name and anything on the right of the
separator will be assumed as the token value. white spaces to the left and right of the
token and the value will be removed by the function. if no separator is found then function
will return token as the whole string with white spaces removed from left and right
token and values are assumed to be long enough to hold the strings extracted
+ added error checks for strings [08/11/03 - Anoop]
+ changed separator to const [08/11/03 - Anoop]
string - the string to be processed
token - the name of the token read from the string
separator - the array of characters to be used as separators.
created 		07-24-03	Anoop John
last Modified	08-11-03	Anoop John
from "HALabUtils.h"
*/
char* getTokenNameAndValue(const char *string, char *token, const char *separator, char *value)
{
	size_t curStrPos;				/*the position in the current string*/

	/*added error check for null string [added 08/11/03 - Anoop]*/
	if(string == 0 || token == 0 || separator == 0 || value == 0)
		errorOut(("Null Pointer passed to function"));
#ifdef VERBOSE
	printf("Extracting token and value from %s\n", string);	
	halt;
#endif	
	curStrPos = 0;
	/*find the separator or end of string*/
	while(string[curStrPos] != 0 && strchr(separator, string[curStrPos]) == 0)
		curStrPos++;
	/*if separator was not found*/
	if(string[curStrPos] == 0)
	{
		strcpy(token, string);
		strcpy(value, "");
	}
	else
	{
		strncpy(token, string, curStrPos);
		/*strncpy doesnt null terminate*/
		token[curStrPos] = 0;
		strcpy(value, &string[curStrPos+1]);
	}
	removeExtraWhiteSpaces(token);
	removeExtraWhiteSpaces(value);
#ifdef VERBOSE
	printf("Token = %s\nValue = %s\n", token, value);
	halt;
#endif
	return token;
}
/*--------------------------------------------------------------------------------*/
/*Function that reads in one data item from a file given the format
 as a string like in the normal scanf function
 fp 		- the pointer to the file from which it should read data from 
 format 	- equivalent to scanf format but only one variable will be read
 target  - address of the variable into which the data will be read into
 created 		:02/01/03	- Anoop John 
 last modified	:02/01/03	- Anoop John
from "HALabUtils.h"
 */
long readDataWithComments(FILE* fp, char*format, void* target)
{
	char tempc, foundeof=0;
	
	/*reads the data*/
	if(fscanf(fp, format, target)!=1)
		errorOut(("Corrupted input file or invalid format string(%s)", format));
	
	/*skips till end of line either CR or LF*/
	tempc=0;
	while((tempc!='\n')&&(tempc!='\r'))
	{
		tempc=fgetc(fp);
		if(foundeof)		/*Added check for ending CR or LF 01-27-03*/
			errorOut(("Missing ending CR or LF character in last line.\nPossible corrupted file "));
		if(feof(fp))
			foundeof=1;
	}
	return 0;
}
/*---------------------------------------------------------------------------------------*/
/* *************************************************************** */
/* CreateFolder(pth, fld_name, mode)       Boyang Li 100909      */
/*---------------------------------------------------------------*/
/* Inputs:                                                       */
/*     pth      -- path of folder (e.g. "/Users/A/Documents/")   */
/*     fld_name -- name of folder (e.g. "untitled folder")       */
/*     mode     -- mode of file (see below)                      */
/*---------------------------------------------------------------*/
/* mode:                                                         */
/* 0 -- S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH (default)                */
/* S_IRWXU -- read, write, execute/search by owner               */
/* S_IRUSR -- read permission, owner                             */
/* S_IWUSR -- write permission, owner                            */
/* S_IXUSR -- execute/search permission, owner                   */
/* S_IRWXG -- read, write, execute/search by group               */
/* S_IRGRP -- read permission, group                             */
/* S_IWGRP -- write permission, group                            */
/* S_IXGRP -- execute/search permission, group                   */
/* S_IRWXO -- read, write, execute/search by others              */
/* S_IROTH -- read permission, others                            */
/* S_IWOTH -- write permission, others                           */
/* S_IXOTH -- execute/search permission, others                  */
/* *************************************************************** */
char * CreateFolder(const char *full_path, int mode_type)
/*
 creates a new folder at desired path
 will check if folder_name exists and if it does will create "foldername_i" where i is a new folder
 no folders will be overwritten
 mode specifies permissions, etc - 0 is all permissions okay
 returns 0 if okay
from "HALabUtils.h"
 */
{
	int status = -1, i = 0, str_len, debug = 2;
	int FolderNameMax = 128;	
	char *pth_f;
	mode_t mode;
	
	if (debug == 1)	{	printf("in function CreateFolder");	halt;	}
	str_len = strlen(full_path) + FolderNameMax;
	pth_f = (char *) memAlloc (str_len, sizeof(char), "CreateFolder: pth_f");	
	sprintf(pth_f, "%s", full_path);
	
	if (mode_type == 0)
		mode = S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH;
	else
	{	printf("CreateFolder() not ready for mode_type != 0");	halt;	}
	
	while (status == -1)
	{
		status = mkdir(pth_f, mode);	
		if (status == -1)
		{
			i++;
			if (i > 20)
			{	printf("all foldernames used through 20 for %s", full_path);	halt;	}
			sprintf(pth_f, "%s_%d", full_path, i);			
		}
		else 
		{
			if (debug > 0)	{	printf("folder: %s created\n", pth_f);	}	
		}
	}
	return pth_f;
}
/*--------------------------------------------------------------------------------*/
/*function to check if folder name ends with the directory separator and if so then
remove the directory separator from the end of the foldername. 
+ added error checks for strings [08/11/03 - Anoop]
folderName - the folderName to be checked
dirSep - the directory separator
created 		07-24-03	Anoop John
last Modified	08-11-03	Anoop John
from "HALabUtils.h"
*/
char* fixFolderName(char *folderName, const char *dirSep)
{
	size_t dirSepLen;				/*length of the directory separator string*/
	size_t folderNameLength;		/*length of the folder name*/

	/*added error check for null string [added 08/11/03 - Anoop]*/
	if(folderName == 0 || dirSep == 0)
		errorOut(("Null Pointer passed to function"));
	folderNameLength = strlen(folderName);
	dirSepLen = strlen(dirSep);
	/*check if the foldername ends with directory separator*/
	if(strstr(&folderName[folderNameLength-dirSepLen], dirSep) != 0)
		folderName[folderNameLength-dirSepLen] = 0;
	return folderName;
}
/*--------------------------------------------------------------------------------*/
char * get_time_date_str(void)
/* uses time.h functions to get current date/time
	str ends with '\n' so this is replaced with '\0'
from "HA0606lib.h"
 */
{
	long eoln_found, linepos;
	char *time_date_str;
	time_t systime; 
	struct tm *currtime;
	systime = time(NULL);
	currtime = localtime(&systime);
	linepos = 0;
	time_date_str = get_substr_str_tag_delim(asctime(currtime), "firstpos", &linepos, "lastpos", 1000);
	
/*	systime = time(NULL);
	currtime = localtime(&systime);
	asctime(currtime);
*/	
	/* replace \n with \0 in time string */
	eoln_found = count_char_in_str(time_date_str, '\n');
	if (eoln_found == 1)
	{
		linepos = get_stringpos(time_date_str, "\n", 0);
		time_date_str[linepos] = '\0';	
	}
	
	return time_date_str;
}
/*--------------------------------------------------------------------------------*/
long get_file_dataline_ct(char readfile[], long *last_line_eoln)
/* 	returns number of lines with alphnum characters in a file.  The last line is counted if it does not end in an '\n' 
	recode for version that skips initial comment lines?
	 from "HA0606lib.h"

*/
{
	long data_line_ct, tot_line_ct, line_char_ct, debug=0;
	FILE * fpreadfile;
	char ch;
	
	fpreadfile = fileOpen(readfile, "r");
 	data_line_ct = tot_line_ct = line_char_ct = 0;
 	while((ch = fgetc(fpreadfile)) != EOF)
 	{
 		if (ch == '\n')
 		{
 			tot_line_ct++;
 			if (line_char_ct > 0)
 				data_line_ct++;
 			line_char_ct = 0;		/* counts number of isalnum char preceding each '\n' */
 		}
 		else
 			if (isalnum(ch) != 0)
 				line_char_ct++;
 	}  
 	/* count the last line if there are alphanum and no '\n' */
 	if (line_char_ct > 0)
 	{
 		data_line_ct++;
 		tot_line_ct++;
 		*last_line_eoln = 0;
/* 		printf("Last line in %s contains alnum char but no line ending char\n", readfile);
 		printf("This will cause problems with get_line_from_file function");
 		printf("Please run file_read_check program");
 		halt; */
 	}
 	else
 	{
/* 		printf("last line in %s with alnum char has a line ending char\n", readfile); */
 		*last_line_eoln = 1;
 	}
 	if (debug == 1)
	{
		printf("%ld lines in %s\n", tot_line_ct, readfile);
		printf("%ld lines with alnum char in %s\n\n", data_line_ct, readfile);
	}
	fclose(fpreadfile);

	return tot_line_ct;
}
/*--------------------------------------------------------------------------------*/
void copy_file_w_time(char readfile[], char copyfile[])
/* 	Reads a file and copies all contents to another file
 The readfile and copy file (including path) are passed to the function. 	
 problem with line endings in copy_files? 
 from "HA0606lib.h"
 */
{
	FILE * fpreadfile, *fpcopyfile;
	long curr_str, line_num, last_line_eoln, eof_found /*, linepos */;
	char *data_line, *time_date_str;
	//	time_t systime; 
	//	struct tm *currtime;
	
	last_line_eoln = -9;
	line_num = get_file_dataline_ct(readfile, &last_line_eoln);	
	
	fpreadfile = fileOpen(readfile, "r");
	fpcopyfile = fileOpen(copyfile, "w");
	
	time_date_str = get_time_date_str();	
	fprintf(fpcopyfile, "/*-  %s  -*/\n", readfile);
	fprintf(fpcopyfile, "/*-  %s  -*/\n", time_date_str);
	for (curr_str = 1; curr_str <= line_num; curr_str++)
	{
		data_line	 	= get_line_from_file_rev(fpreadfile, readfile, &eof_found);
		fprintf(fpcopyfile, "%s", data_line);
		if ((curr_str < line_num) || (last_line_eoln == 1))
			fprintf(fpcopyfile, "\n");
		free(data_line);
	}
	fclose(fpreadfile);
	fclose(fpcopyfile);
	
	return;
}
/*--------------------------------------------------------------------------------*/
char *get_line_from_file_rev( FILE *fptr, char *filename, long *eof_found)
/* Function to get line from a file. It accepts file pointer and file name as input. It counts the number of characters
 in the line excluding the '\n' character. It then allocates count+1 number of memory spaces to a string and reads in the 
 line character by character into the string. It terminates the string with the '\0' character. It returns the next line 
 of the file as string . This function does not, however, free the memory after allocating it.
 if the next line is the end of file, the function returns 0 
 the "rev" version of this function reads the last line of a file that does not include a '\n' (line ends with EOF)	
 note that empty lines return a string of length 1 with string[0] = '\0' 
  from "HA0606lib.h"
*/
{
   char *strng, ch;
   long count = 0, i, go_back_char_ct, debug = 0;
      
	/* Count the number of characters into the variable "count", excluding the '\n' character so that, that many memory spaces could be allocated for the string */       
 	*eof_found = 0;
 	while(((ch = fgetc(fptr)) != '\n') && (*eof_found == 0))
 	{
		if (ch == EOF) /* If no new lines are found till end of file, report failure */
			*eof_found = 1;
 		else
			count++;
	}      
      	
	/* Return file pointer to starting of the line by moving (count + 1) times backwards (one is added to include the '\n' character but not EOF) */     
	go_back_char_ct = count;
	if (*eof_found == 0) 
		go_back_char_ct++;
  	if ((fseek(fptr, -((go_back_char_ct)*sizeof(char)), SEEK_CUR)) != 0)
    	{	printf("Error seeking the position in get_line_from_file_rev for file %s", filename);	halt;	}
    
	/* Allocate memory for the string equal to (count+1) times the size of a character to include '\0' character at the end of the string */      
    strng = (char *) memAlloc ((count + 1), sizeof(char), "get_line_from_file: strng");

	/* Read the line character by character */     
    for(i = 0; i <= (count - 1); i++)
  		if ((fscanf(fptr, "%c", &strng[i])) != 1)	/* Error checking to see if each character has been read correctly */
  			{	printf("Error reading the character %c in the file %s", strng[i], filename);	halt; }
 
    if ((ch = fgetc(fptr)) != '\n') /* To make the file pointer point to the character after the '\n' */  
    	if (ch != EOF)	
			{	printf("\n File pointer error in %s\nProblem with the end of the line.", filename);	halt; }

    strng[count]='\0';
	
	if (debug == 1)
		if (*eof_found == 1)
		{
			printf("\t\t\tget_line_from_file_rev(): eof_found in -%s-\nstring: -%s-\n", filename, strng);
			halt;
		}
    
    return(strng);   
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------*/
long get_stringpos(char *parent, char *tag, long start_index)
/*  returns the index of the first position of tag string within parent following (and including the start_index)
	returns -9 if the tag is not found 
	from "HA0606lib.h"
*/
{
	long pos, i;
	char *parent_mod, *tag_ptr;

	if (start_index >= strlen(parent))	/* start_index is location of string termination or beyond */
		return -9;
	
	parent_mod = parent; 				/* parent_mod will be a pointer to first char in parent to search from */
	for (i=1; i<=start_index; i++)
		parent_mod++;
		
	tag_ptr = strstr(parent_mod, tag); 	/* tag_ptr points to first location of tag within the parent string after position=start_index */
	if (tag_ptr == NULL)
		return -9;
	pos = tag_ptr - parent; 			/* pos is a long, tag_ptr and parent are pointers to char, so pos is the distance from the beginning of parent to the beginning of tag_ptr */
	
/*printf("parent: \n%s\ntag: %s\nstart_index: %5ld\npos: %5ld\n", parent, tag, start_index, pos);
*/	
	return pos;
}	
/*---------------------------------------------------------------------------------------*/
long count_char_in_str(char *input_string, char char_to_count)
/*
 from "HA0606lib.h"
*/
{
	long count, pos_ct, i;
	
	if (input_string == NULL)
		return -9;
	count=0;
	pos_ct = strlen(input_string) - 1;  /* pos_ct is the last index before the string terminating '/0' */
	for (i=0; i<=pos_ct; i++)
		if (input_string[i] == char_to_count)
			count++;

	return count;
}
/*----------------------------------------------------------------------------------------------------*/
char * get_substr_str_tag_delim(char *parent_str_orig, char *tag, long *start_index, char *delimiter_str, long max_substring)
/*  returns a pointer to string for char found between "tag" and "delimiter" strings within a parent string
	the new string is malloc-ed
	start_index determines where in the parent string to begin search
		it is passed as the address of a long value
		the value is changed to point to the index of the last position of the delimiter (this is important! the same delimiter can be used as the tag for the next target)
		if tag = "firstpos", then use index start_index for beginning of target_string rather than search for tag_string
		if delimiter = "lastpos", then use last position in string as end of target string rather than search for delimiter_string
	if the found string is >max_substring char, a warning is printed and the excutable halts (this can be used to check for errors in the readfile)
	function gives a warning, halts, and returns NULL pointer if tag is not found in parent
 from "HA0606lib.h"
*/
{
	char *target_string, *parent_str;
	long tag_startpos, tag_len, target_len, target_startpos, target_endpos, i, j;
	
	parent_str = parent_str_orig;			/* parent string points to position "index" in the original parent string */
	for (i=1; i<= *start_index; i++)
		parent_str++;
	
	if (strcmp(tag, "firstpos") != 0)
	{
		tag_startpos = get_stringpos(parent_str, tag, 0);	/* tag_startpos is the index of the first char following the tag string in parent_str */
		if (tag_startpos == -9)	/* tag not found */
		{	
			printf("get_substring_tag_delim: -%s- not found in -%s-\nreturning NULL pointer", tag, parent_str); 
			halt;	
			return NULL;	
		}

		tag_len = strlen(tag); 	/* tag_len is the length of the string tag excluding the '\0' character */
		
		target_startpos = tag_startpos;	/* increment target_startpos to first position after the tag string */
		for (i=1; i<= tag_len; i++)
			target_startpos++;
	}
	else
		target_startpos = 0;

	if (strcmp(delimiter_str, "lastpos") != 0)
	{
		target_endpos = get_stringpos(parent_str, delimiter_str, target_startpos); /* get last index before delimiter */
		if (target_endpos == -9)
			{printf("get_substring_tag_delim: -%s- not found in -%s- after pos %ld\n", delimiter_str, parent_str, target_startpos); halt;}
		target_endpos--;	
	}
	else 
		target_endpos = strlen(parent_str) - 1;
		
	target_len = target_endpos - target_startpos + 1;	/* number of char that will be in target_str - the value will be 0 if there is are no char in the target */
	if (target_len > max_substring)
		{	printf("delimiter -%s- not found within %0ld of tag: -%s- for -%s-\n", delimiter_str, max_substring, tag, parent_str);	halt;	}
	
	target_string = (char *) memAlloc ((target_len + 1), sizeof(char), "get_substring_tag_delim: target_string"); 
    if (target_string == NULL)
    	{	printf("\nMemory allocation failed for parent string %s and tag %s", parent_str, tag);	halt;	}
	 		 
	for(j=0, i=target_startpos; i <= target_endpos; i++, j++)	/* fill target_string */
      	target_string[j] = parent_str[i];
	target_string[target_len] = '\0'; 			/* if target_len_0, then the termination will be the first char */
	
	if (strlen(target_string) != target_len)
	{	printf("\nget_substring_tag_delimstrlen(target_string) != target_len.  May be bug in code.\nparent: %s\ntarget: %s", parent_str, tag);	halt;	}
	
	if (strcmp(delimiter_str, "lastpos") != 0)	/* start_index value becomes the last index of the end delimiter in the original parent string */
		*start_index = *start_index + target_endpos + strlen(delimiter_str);	/* start_index value is now the index of the last position of the delimiter */
	else
		*start_index = *start_index + target_endpos + 1; /* should point to '\0' at end of string */

	return(target_string);
}
/*****************************************************************************************/
long StringToLong(char * InputString)
/*
from "IO_functions.h"
*/
{
	long pos=0;
	long Length;
	long Number;

	while(InputString[0] == ' ' || InputString[0] == '\t' || InputString[0] == '\n')
		InputString++;	
	Length = (long) strlen(InputString);
	if(Length == 0)
		errorOut(("Null string or a non-valid string passed to convert to long"));
	else if (Length >= 10)
	{
		printf("The number read was %s. Please quit if the number is not supposed to be this large\n", InputString);
		halt;
	}	
	while(InputString[pos] != '\0' && InputString[pos] != '\n' && pos<Length)
	{
		if(isspace(InputString[pos]))
		{
			pos++;
			continue;
		}
		if ((isdigit(InputString[pos]) == 0) && (InputString[pos] != '-'))	/* use isdigit instead of isalnum returns 8 when the character is a numeric character */
			errorOut(("Could not convert number to long v2\nRead: $%s$\n", InputString));
		pos++;
	}
	Number = (long) atoi(InputString);
	if(errno == ERANGE)
		errorOut(("Error in converting the string %s to long\n",InputString));
	return Number;
}
/*****************************************************************************************/
double StringToDouble(char * InputString)
/*
from "IO_functions.h"
*/
{
	long pos=0;
	long Length;
	double Number;
	
	Length = strlen(InputString);
	if(Length == 0)
		errorOut(("Null string or a non-valid string passed to convert to double"));
	while(InputString[pos] != '\0' && InputString[pos] != '\n' && pos<Length)
	{
		if(isspace(InputString[pos]))
		{
			pos++;
			continue;
		}
		if(isdigit(InputString[pos]) == 0 && InputString[pos] != '-' && InputString[pos] != '+' && InputString[pos] != '.' && InputString[pos] != 'e' && InputString[pos] != 'E')	
			/* isalnum returns 8 when the character is a numeric character */
			errorOut(("Could not properly convert the number to double\nFunction received: $%s$\n", InputString));
		pos++;
	}
	Number = (double) atof(InputString);
	if(errno == ERANGE)
		errorOut(("Error in converting the string %s to double\n",InputString));
	return Number;
}/*---------------------------------------------------------------------------------------*/
unsigned long programTimer(long m)
/* function to find the time taken inside program.
m = 0 		initialize timer
m = 1 		print out time since initialization
m = 2		print out time since last call
m = 3		return time since initialization
m = 4		return time since last call
Created 		06/06/03	Anoop John
Last Modified	06/17/03	Anoop John 
from "HALabUtils.h"
*/
{
	static time_t start=0;			/*the time when the timer was initialized*/
	static time_t end=0; 			/*the time when the function was last called*/
	time_t sttemp=0;
	time_t diff;
	if(start == 0 && m != 0)
		errorOut(("Timer not initialized. initialize with programTimer(0)"));
	switch(m)
	{
		/*initialize timer*/
		case 0:
		{
			start = time(NULL);
			break;
		}
		/*print time from start*/
		case 1:
		{
			end = time(NULL);
			diff = difftime(end, start);
			printf("\nTotal Time from start:  %02ld:%02ld:%02ld = %02ld:%02ld\n", diff/3600, (diff%3600)/60, diff%60, diff/60, diff%60);
			break;
		}
		/*print time since last call*/
		case 2:
		{
			if(end)
				sttemp = end;
			else
				sttemp = start;
			end = time(NULL);
			diff = difftime(end, sttemp);
			printf("Last Interval:  %02ld:%02ld:%02ld = %02ld:%02ld\n", diff/3600, (diff%3600)/60, diff%60, diff/60, diff%60);
			break;
		}
		/*return time in sec from start*/
		case 3:
		{
			end = time(NULL);
			diff = difftime(end, start);
			return diff;
			break;
		}
		/*return time in sec since last call*/
		case 4:
		{
			if(end)
				sttemp = end;
			else
				sttemp = start;
			end = time(NULL);
			diff = difftime(end, sttemp);
			return diff;
			break;
		}
		default:
		{
			errorOut(("Invalid parameter passed to programTimer. (Valid range = 0-4)"));
			break;
		}
	}
	return 0;
}
/*---------------------------------------------------------------------------------------*/
/*function to handle memory allocation, reallocation  and freeing so that memory errors are
avoided  and thus  preventing silent  crashes. function keeps an internal array of pointers
and  an array  of  memory allocated  to each pointer.  freeing or reallocating pointers not
allocated through the function will result in error. All memory operatins should ideally be
handled by this function. if a pointer allocated through a normal call to malloc, calloc or
realloc  is passed to this  function then function  will give error. if a pointer allocated
through the function is freed using a  call to free then the internal list for the function 
will no longer be correct and this can lead to errors. memory allocated through this should
always be freed using this function
+ added keeping track of variable names so that in the final printMemoryUsageStatistics the 
	outstanding pointers can actually be located [08/21/03 - Anoop]
+ added option to print summary alone [08/27/03 - Anoop]
+ added peakMemoryUsage output with the summary [09/09/03 - Anoop]
memSize 		- the size of the block to be allocated or reallocated
orig			- the original pointer. for allocateMemory request this will be ignored
variableName	- the name of the variable passed to the function, will be used to keep track of 
				  variables for which memory was allocated
request			- has 6 options currently
				1 = allocateMemory
				2 = reAllocateMemory
				3 = freeMemory
				4 = printStatistics
				5 = printSummary	
				6 = resetFunction	
				7 = freeAllMemory
created 		08-14-03	Anoop John
last Modified	09-09-03	Anoop John
*/
void* memoryHandler(size_t memSize, void *orig, const char *variableName, char request)
{
	static unsigned long numBlocks = 0;
	static void **blockPointers = 0;
	static unsigned long *blockSizes = 0;
	static char **blockNames = 0;
	static unsigned long totalAllocateCalls = 0;
	static unsigned long totalReAllocateCalls = 0;
	static unsigned long totalFreeCalls = 0;
	static unsigned long memoryOutstanding = 0;
	static unsigned long totalAllocatedMemory = 0;
	static unsigned long peakMemoryUsage = 0;
	static unsigned long totalFreedMemory = 0;
	static unsigned long curMaxNumBlocks = 0;
	const unsigned long incrementBlockNum = 1000;
	const char allocateMemory = 1;
	const char reAllocateMemory = 2;
	const char freeMemory = 3;
//	const char printStatistics = 4;
//	const char printSummary = 5;
//	const char resetFunction = 6;
//	const char freeAllMemory = 7;
	unsigned long curBlock;
	int maxStrLen;
	unsigned long curStrLen;
	unsigned long curChar;
	void *newBlock;
	
	/*if request is to reallocate memory and the memSize is 0 then it is same as free call*/
	if(request == reAllocateMemory && memSize == 0)
		request = freeMemory;
	/*if request is to reallocate memory and the original pointer is 0 then it is same as allocateMemory*/
	if(request == reAllocateMemory && orig == 0)
		request = allocateMemory;
	
	switch (request)
	{	
		/*allocate memory*/
		case 1:
				totalAllocateCalls++;
				/*allocate memory for the request*/
				newBlock = calloc(1, memSize);
				if(newBlock != 0)
				{
					/*increment the number of blocks allocated*/
					numBlocks++;
					/*if the number of blocks is greater than the current maximum then allocate more memory*/
					if(numBlocks > curMaxNumBlocks)
					{
						curMaxNumBlocks += incrementBlockNum;
						blockPointers = realloc(blockPointers, sizeof(void*)*(curMaxNumBlocks+1));
						blockSizes = realloc(blockSizes, sizeof(unsigned long)*(curMaxNumBlocks+1));
						blockNames = (char **)realloc(blockNames, sizeof(char*)*(curMaxNumBlocks+1)); 
						if(blockPointers == 0 || blockSizes == 0 || blockNames == 0)
							errorOut(("Memory Allocation error. Couldnt reallocate memory for the internal arrays"));
					}
					/*add pointer and size to list*/
					blockPointers[numBlocks] = newBlock;
					blockSizes[numBlocks] = memSize;
					blockNames[numBlocks] = calloc(strlen(variableName)+1, sizeof(char));
					if(blockNames[numBlocks] == 0)
						errorOut(("Memory Allocation error. Couldnt reallocate memory for the variable name"));
					strcpy(blockNames[numBlocks], variableName);
					memoryOutstanding += memSize;
					totalAllocatedMemory += memSize;
					if(peakMemoryUsage < memoryOutstanding)
						peakMemoryUsage = memoryOutstanding;
				}
				return newBlock;
				break;
		/*if request is to reallocate memory*/
		case 2:
				totalReAllocateCalls++;
				/*check if the pointer is already in the list*/
				for(curBlock = 1; curBlock <= numBlocks; curBlock++)
					if(blockPointers[curBlock] == orig)
						break;
				/*if the pointer is not in the list then error*/
				if(curBlock > numBlocks)
					errorOut(("Cannot reallocate memory for pointer(%p) passed to function. The pointer was not allocated through memRealloc or memAlloc", orig));
				/*try to reallocate memory if it is already in the list*/
				newBlock = realloc(orig, memSize);
				/*if successfully reallocated memory then add pointer and size to list*/
				if(newBlock != 0)
				{
					/*add pointer and size to list*/
					blockPointers[curBlock] = newBlock;
					memoryOutstanding -= blockSizes[curBlock];
					totalFreedMemory += blockSizes[curBlock];
					free(blockNames[curBlock]);
					blockNames[curBlock] = calloc(strlen(variableName)+1, sizeof(char));
					if(blockNames[curBlock] == 0)
						errorOut(("Memory Allocation error. Couldnt reallocate memory for the variable name"));
					strcpy(blockNames[curBlock], variableName);
					memoryOutstanding += memSize;
					blockSizes[curBlock] = memSize;
					totalAllocatedMemory += memSize;
					if(peakMemoryUsage < memoryOutstanding)
						peakMemoryUsage = memoryOutstanding;
				}
				return newBlock;
				break;
		/*if request is to free memory*/
		case 3:
				totalFreeCalls++;
				/*check if the pointer is already in the list*/
				for(curBlock = 1; curBlock <= numBlocks; curBlock++)
					if(blockPointers[curBlock] == orig)
						break;
				/*if the pointer is not in the list then error*/
				if(curBlock > numBlocks)
					errorOut(("Cannot free memory for pointer(%p) passed to function. The pointer was not allocated through memRealloc or memAlloc", orig));
				/*if already in the list free memory and update the list*/
				free(orig);
				blockPointers[curBlock] = blockPointers[numBlocks];
				memoryOutstanding -= blockSizes[curBlock];
				totalFreedMemory += blockSizes[curBlock];
				blockSizes[curBlock] = blockSizes[numBlocks];
				free(blockNames[curBlock]);
				blockNames[curBlock] = blockNames[numBlocks];
				numBlocks--;
				return 0;
				break;
		/*if request is to print statistics*/
		case 4:
				/*get length of strings*/
				maxStrLen = 0;
				for(curBlock = 1; curBlock <= numBlocks; curBlock++)
				{
					curStrLen = strlen(blockNames[curBlock]);
					if(maxStrLen < curStrLen)
						maxStrLen = curStrLen;
				}
				maxStrLen++;
				printf("=====================================\n");
				printf("===memoryHandler Statistics==========\n");
				if(numBlocks > 0)
				{
					printf("=====================================\n");
					printf("%-18s:%10s:%*s\n", "Outstanding block", "Size", maxStrLen, "Name");
					printf("------------------:----------:");
					for(curChar = 1; curChar <= maxStrLen; curChar++)
						printf("-");
					printf("\n");
					for(curBlock = 1; curBlock <= numBlocks; curBlock++)
						printf("%-18p:%10lu:%*s\n", blockPointers[curBlock], blockSizes[curBlock], maxStrLen, blockNames[curBlock]);
					printf("-----------------------------:");
					for(curChar = 1; curChar <= maxStrLen; curChar++)
						printf("-");
					printf("\n");
				}
				/*dont break here continue on to print summary*/
		/*print the summary*/
		case 5:
				printf("=====================================\n");
				printf("Summary\n");
				printf("=======\n");
				printf("Total pointers outstanding : %lu\n", numBlocks);
				printf("Total memory outstanding   : %lu\n", memoryOutstanding);
				printf("Total memory allocated     : %lu\n", totalAllocatedMemory);
				printf("Total memory freed         : %lu\n", totalFreedMemory);
				printf("Number of allocate calls   : %lu\n", totalAllocateCalls);
				printf("Number of reAllocate calls : %lu\n", totalReAllocateCalls);
				printf("Number of free calls       : %lu\n", totalFreeCalls);
				printf("Peak memory usage          : %lu\n", peakMemoryUsage); 
				printf("=====================================\n");
				return 0;
				break;
		/*clean memory allocated locally in function*/
		case 6:
				if(numBlocks != 0)
				{
					printf("You are about to clear the internal lists in the function memoryHandler before all the outstanding memory is freed. Free memory allocated using calls to free()");
					halt;
				}
				if(blockPointers)
					free(blockPointers);
				if(blockSizes)
					free(blockSizes);
				if(blockNames)
				{
					for(curBlock = 1; curBlock <= numBlocks; curBlock++)
						if(blockNames[curBlock])
							free(blockNames[curBlock]);
						else
							errorOut(("blockNames[curBlock] cannot be zero"));
					free(blockNames);
				}
				blockPointers = 0;
				blockSizes = 0;
				blockNames = 0;
				numBlocks = 0;
				return 0;
				break;
		/*free all memory allocated through this function*/
		case 7:
				for(curBlock = 1; curBlock <= numBlocks; curBlock++)
				{
					if(blockPointers[curBlock])
						free(blockPointers[curBlock]);
					else
						errorOut(("blockPointers[curBlock] cannot be 0"));
					if(blockNames[curBlock])
						free(blockNames[curBlock]);
					else
						errorOut(("blockNames[curBlock] cannot be zero"));
				}
				free(blockNames);
				free(blockPointers);
				free(blockSizes);
				blockPointers = 0;
				blockSizes = 0;
				numBlocks = 0;
		/*if the request passed is not one identified by the function error*/
		default:
				errorOut(("Request has to be one of allocateMemory(1), reAllocateMemory(2), freeMemory(3), printStatistics(4), resetFuncion(5) or freeAllMemory(6)"));
				break;
	}
	return 0;
}
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*Multinomial deviate function. Function fills in the array nn with the number of 
 successes for the corresponding probability of successes given in the inn array and
 the total number of experiments n. k is the size of the two arrays.(1 based arrays)
 Function Created:  Hiroshi Akashi
 */
void multdev(double inn[], long k, long n, long nn[])
{
	long i, r;
	double sum;
	
	sum = 0.0;
	for(i=1;i<=k;i++)
	{
		if( (n > 0) && (sum < 1.0) )
		{
			r = (long)(bnldev(inn[i]/(1.0 - sum), n));
			nn[i] = r;
			sum += inn[i];
			n -= r;
		}
		else
		{
			nn[i] = 0;
		}
	}
	return;
}
/*--------------------------------------------------------------------------------*/
/*Function returns the binomial deviate. Given the probability pp and the number of
 experiments the function returns the number of successes(long) as a double value.
 Function taken from Numerical Recipes for C IIed.
 */
double bnldev(double pp, long n)
{
	long j;
	static long nold=(-1);
	double am, em, g, angle, p, bnl, sq, t, y;
	static double pold=(-1.0), pc, plog, pclog, en, oldg;
	
	p=(pp <= 0.5 ? pp : 1.0-pp);
	am=n*p;
	if (n < 25) {
		bnl=0.0;
		for (j=1;j<=n;j++)
			if (ran2() < p) ++bnl;
	} else if (am < 1.0) {
		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= ran2();
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {
		if (n != nold) {
			en=n;
			oldg=gammln(en+1.0);
			nold=n;
		} if (p != pold) {
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc);
		do {
			do {
				angle=PI*ran2();
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0));
			em=(long)(em);
			t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
										  -gammln(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (ran2() > t);
		bnl=em;
	}
	if (p != pp) bnl=n-bnl;
	return bnl;
}
/*--------------------------------------------------------------------------------*/
/*Modified Ran2 function. Original function from the Numerical Recipes in C IIed. 
 This function uses a global idum. returns random number in the range (0, 1].
 */
#define IM1 2147483563
#define IM2 2147483399
#define AM2 (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 256
#define NDIV2 (1+IMM1/NTAB)
/*
 uncomment these lines as well as the ones at the end
 of the function to get (0, 1) output.
 */
/*
 #define EPS DBL_EPSILON
 #define RNMX (1.0-EPS)
 */

double ran2(void)
{
	long j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	
	if (*idum <= 0) 
	{ 
		if (-(*idum) < 1) *idum=1; /*Be sure to prevent idum = 0.*/
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+32;j>=0;j--) 
		{
			/*Load the shuffle table (after 32 warm-ups).*/
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1; /*Start here when not initializing.*/
	/*Compute idum=(IA1*idum) % IM1 without*/
	/*overflows by Schrage's method. */
	*idum=IA1*(*idum-k*IQ1)-k*IR1; 
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	/*Compute idum2=(IA2*idum) % IM2 likewise.*/
	idum2=IA2*(idum2-k*IQ2)-k*IR2; 
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV2; /*Will be in the range 0..NTAB-1.*/
	iy=iv[j]-idum2; 
	/*Here idum is shuffled, idum and idum2 are
	 combined to generate output. */
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	/*the output will be (0, 1]*/
	return  (double)AM2*iy;
	/*	
	 Here the output will be (0, 1)
	 if ((temp=AM2*iy)>RNMX) return RNMX
	 else return temp;
	 */
	
}
#undef IM1
#undef IM2
#undef AM2
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV2
/*--------------------------------------------------------------------------------*/
/*Function returns the poisson deviate given a mean value xm. Original function from 
 the Numerical Recipes in C IIed. modified to use global idum 
 */
double poidev(double xm)
{
	static double sq, alxm, g, oldm=(-1.0);
	double em, t, y;
	
	if (xm < 12.0) 
	{
		if (xm != oldm) 
		{
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do 
		{
        	++em;
        	t *= ran2();
		} while (t > g);
	} 
	else 
	{
		if (xm != oldm) 
		{
        	oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-gammln(xm+1.0);
		}
		do 
		{
        	do 
        	{
				y=tan(PI*ran2());
				em=sq*y+xm;
			} while (em < 0.0);
			em=(long)(em);
    		t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (ran2() > t);
	}
	return em;
}
/*--------------------------------------------------------------------------------*/
/*ln gamma(x) function from the Numerical Recipes in C IIed. */
double gammln(double xx)
{
	double x, tmp, ser;
	static double cof[6]={76.18009173, -86.50532033, 24.01409822, -1.231739516, 0.120858003e-2, -0.536382e-5};
	long j;
	
	x = xx-1.0;
	tmp = x + 5.5;
	tmp -= (x+0.5) * log(tmp);
	ser = 1.0;
	for (j=0; j<=5; j++) 
	{
		x += 1.0;
		ser += cof[j]/x;
	}
	return -tmp+log(2.50662827465*ser);
}
/*--------------------------------------------------------------------------------*/
/*function to find the time. starts if called with 
 parameter 0. if called with parameter 1 it gets the 
 total time and exits. if called with parameter 2 it will 
 find the time since the last call and keep the counter on*/
long getTime(long m)
{
	static long start=0, end=0;
	static long sttemp=0;
	long diff;
	
	if(m==0)
	{
		start=time(NULL);
	}
	/*time from start*/
	else if(m==1)
	{
		end=time(NULL);
		diff=difftime(end, start);
		printf("\nTotal Time from start:  %02ld:%02ld:%02ld = %02ld:%02ld\n", diff/3600, (diff%3600)/60, diff%60, diff/60, diff%60);
	}
	/*time since last call*/
	else if(m==2)
	{
		if(end)sttemp=end;
		else sttemp=start;
		end=time(NULL);
		diff=difftime(end, sttemp);
		printf("Last Interval:  %02ld:%02ld:%02ld = %02ld:%02ld\n", diff/3600, (diff%3600)/60, diff%60, diff/60, diff%60);
	}
	/*time from start*/
	else if(m==3)
	{
		end=time(NULL);
		diff=difftime(end, start);
		return diff;
	}
	/*time since last call*/
	else if(m==4)
	{
		if(end)sttemp=end;
		else sttemp=start;
		end=time(NULL);
		diff=difftime(end, sttemp);
		return diff;
	}
	return 0;
}
/*---------------------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------------*/
/*Function to compare to long values given that they are passed
as constant void pointers. Just subtracts and returns values. To be used
with qsort*/
int  compareUL(const void*a, const void*b)
{
	if((*(long*)a)>(*(long*)b))
		return 1;
	else if((*(long*)a)<(*(long*)b))
		return -1;
	return 0;
}
/*---------------------------------------------------------------------------------------*/


