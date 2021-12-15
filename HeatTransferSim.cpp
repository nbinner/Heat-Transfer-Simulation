/************************************************************************************************************
Mech 7171 Engineering Programing, Fall 2021
Capstone Project - 2D Conductive Heat Transfer on a Uniform Flat Plate

Purpose: The purpose of this lab is to solve the heat conduction equation in 2D using Finite-Difference
		 (F-D) methods.  The solution will be compared to an analytic solution and solution
		 convergence will be tracked.  This lab will exercise all facets of C taught in Mech 7171.

Author(s):     Muneer Almasyabi, Nathan Binner, Alex Sung
Student ID(s): A01061394, A01159743, A01163512

Declaration:
We, Munner, Nathan and Alex declare that the following program was written by me/us.

Date Created: 2021-11-22
************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

//------- GLOBAL CONSTANTS ----------------------------------------------------------------------------------
const double T0 = 0.0;     // normal background wall temperature (for initializing!)

const char* SIMULATION_DATA_FILE_HEADER1 = "Simulation    w      h      dx      dy";
const char* SIMULATION_DATA_FILE_HEADER2 = "Boundary Conditions";
const char* ENDMESSAGE = "Thank you for running this program!";
const char* DASHES = "--------";
const char* SIMULATIONS_INPUT_DATA_FILE = "simulations.in";
const char* SIMULATION_CHOICE = "Select Simulation Case:";

const int CASE_TYPE_A = 0;
const int CASE_TYPE_B = 1;
const int CASE_TYPE_C = 2;
const int CASE_TYPE_TEST = 3;
const int MAX_CASE_NAME_SIZE = 40;

const int TABLE_COLUMN_WIDTH = 50;
const int TABLE_MARGIN_SIZE = 3;

const int NUM_WALLS = 4;
const int TOP = 0;
const int BOTTOM = 1;
const int LEFT = 2;
const int RIGHT = 3;

const int BC_TYPE_CONST = 0;
const int BC_TYPE_COSINE = 1;
const int BC_TYPE_INSULATED = 2;
const int BC_TYPE_POLY = 3;
const int BC_TYPE_SINE = 4;

const double PI = 3.141592653589793;
const int MAX_BUFF_SIZE = 1024;              // for reading lines from a file
const int MAX_ITER = 1000000;                // maximum iterations for F-D
const double MAX_RESIDUAL = 1.0e-8;          // solution accuracy
const double MAX_TRANSIENT_RESIDUAL = 0.01;  // transient solution accuracy

//--- Table Border Characters
const unsigned char HL = 196;  // horizontal border line
const unsigned char VL = 179;  // vertical border line
const unsigned char TL = 218;  // top left border symbol
const unsigned char TC = 194;  // top center border symbol
const unsigned char TR = 191;  // top right border symbol
const unsigned char CL = 195;  // left center border symbol
const unsigned char CC = 197;  // center center border symbol (cross)
const unsigned char CR = 180;  // right center border symbol
const unsigned char BL = 192;  // bottom left border symbol
const unsigned char BC = 193;  // bottom center border symbol
const unsigned char BR = 217;  // bottom right border symbol


//------- STRUCTURE DEFINITIONS -----------------------------------------------------------------------------
typedef struct PLATEPOINT
{
	double x, y;   // grid node physical position on the plate
	double T_a;    // grid node temperature (analytic solution)
	double T_fd;   // grid node temperature (finite-difference solution)
	double res;    // grid node residual for finite-difference solution
}
PLATEPOINT;

typedef struct BOUNDARY_CONDITION_DATA
{
	int    nType;     // CONST, SINE, COSINE, POLY, INSULATED
	double Ta, Tb;    // if only one then Ta is used
	double za, zb;    // range (z is one of x or y)
	double ma, mb;    // for BC_TYPE_POLY only
	double k;         // for BC_TYPE_SINE only
}
BOUNDARY_CONDITION_DATA;

typedef struct SIMULATION_DATA    // holds data for each simulation
{
	double w, h, dx, dy;                   // plate width, plate height; x, y cellSizes
	size_t I, J;                           // number of nodes in x and y directions
	int nCaseType;                         // for chosing boundary conditions
	char strCase[MAX_CASE_NAME_SIZE];      // case name
	BOUNDARY_CONDITION_DATA bc[NUM_WALLS]; // one for each wall
}
SIMULATION_DATA;


//------------------------- FUNCTION PROTOTYPES -------------------------------------------------------------
int  nint(double);                           // get the nearest integer to a double value
void waitForEnterKey();                      // robust version of getchar (cleans input buffer)
bool flushInputBuffer2();                    // return true if has garbage left in input buffer
void endProgram(const char*);               // terminates program with optional message
void printRepeatedChar(unsigned char, int);  // prints a character repeatedly
void removeNewline(char* str);               // removes newline at end of string
bool isBlankLine(const char*);              // checks if a line contains only whitespace chars
SIMULATION_DATA* GetSimulationData(SIMULATION_DATA*, int*); // reads a input file to obtain simulation data
int caseTypetoInt(char*);                                   // converts string caseType to an integer
int getUserSimulationChoice(SIMULATION_DATA*, int);         // gets the users sim choice for processing 
int  printHorizontalBorder(char, char);  // Prints the top or bottom border of the array display box
void drawStringLine(const char*, int); // draws each string line within the menu block
void GetCaseAAnalyticalSolution(PLATEPOINT**, const SIMULATION_DATA*); // xmas present! Thanks Dave!
void GetCaseBAnalyticalSolution(PLATEPOINT**, const SIMULATION_DATA*); // xmas present! You're a cool dude
void GetCaseCAnalyticalSolution(PLATEPOINT**, const SIMULATION_DATA*); // xmas present! Appreciate it 
void GetNumericalSolution(PLATEPOINT**, const SIMULATION_DATA);  // numerically calculates the solution of each case
void printSolution(PLATEPOINT**, const SIMULATION_DATA*); // 2nd xmas present!  Prints contour plot data.
PLATEPOINT** initialize(int, SIMULATION_DATA*, PLATEPOINT**); // initializes the dynamic array
PLATEPOINT** SetBoundaryConditions(PLATEPOINT**, SIMULATION_DATA*, int); // sets boundary conditions for each wall
void FreeMemory(PLATEPOINT**, SIMULATION_DATA*, size_t); // frees the memory of the dynamically allocated arrays


//-----------------------------------------------------------------------------------------------------------
int main()
{
	int iS = -1, NS = -1;         // chosen simulation index, number of simulations
	PLATEPOINT** P = NULL;        // For 2D dynamically allocated array for a simulation
	SIMULATION_DATA* SD = NULL;   // the array to hold simulation data for all cases in simulations.in

	SD = GetSimulationData(SD, &NS);
	iS = getUserSimulationChoice(SD, NS);
	P = initialize(iS, SD, P);
	P = SetBoundaryConditions(P, SD, iS);
	GetNumericalSolution(P, SD[iS]);
	if (strcmp(SD[iS].strCase, "A-1") == 0 || strcmp(SD[iS].strCase, "A-2") == 0)
		GetCaseAAnalyticalSolution(P, &SD[iS]);
	else if (strcmp(SD[iS].strCase, "B-1") == 0 || strcmp(SD[iS].strCase, "B-2") == 0)
		GetCaseBAnalyticalSolution(P, &SD[iS]);
	else if (strcmp(SD[iS].strCase, "C-1") == 0 || strcmp(SD[iS].strCase, "C-2") == 0 || strcmp(SD[iS].strCase, "C-3") == 0)
		GetCaseCAnalyticalSolution(P, &SD[iS]);
	printSolution(P, &SD[iS]);
	FreeMemory(P, SD, SD[iS].I);
	waitForEnterKey();
	
	endProgram(NULL);
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  This function reads from the input file to store the data uesd for the simulation
// ARGUMENTS:    *SD: the SIMULATION_DATA array, 
// RETURN VALUE: SD dynamic array, int NS
//               
//               
SIMULATION_DATA* GetSimulationData(SIMULATION_DATA* SD, int* NS)
{
	FILE* fin;                         // file input stream
	errno_t err;                       // error stream
	const char* caseSeps = " \t\n/";   // word delimiters
	char data[MAX_BUFF_SIZE] = {};     // temporary storage array for the tokenized string
	char* tok = NULL, * nextToken = NULL, * pRet = NULL, * pGarbage; // tok variables to tokenize
	int n, i;                          // counters
	int N = 0; // int variable for NS
	// opens the input file and checks to see if the file exists
	err = fopen_s(&fin, SIMULATIONS_INPUT_DATA_FILE, "r"); 
	if (err != 0 || fin == NULL) // checks to see if file exists
	{
		printf("Cannot open file :/"); 
		waitForEnterKey();
		exit(EXIT_FAILURE);
	}
	// reads the input file per line until it reaches the sim data header or the dashes lines 
	while (fgets(data, MAX_BUFF_SIZE, fin) != NULL) 
	{
		// if the current line read is the sim data header, it will skip past it and not do anything
		if (strstr(data, SIMULATION_DATA_FILE_HEADER1) != NULL) continue; 
		// if the current line read is the dashed lines, it will break out the loop 
		if (strstr(data, DASHES) != NULL)break;
	}
	// continues to read through the document line by line
	while (fgets(data, MAX_BUFF_SIZE, fin) != NULL) 
	{
		// checks to see if the current line is a blank line and if it is it will skip it
		if (isBlankLine(data)) continue;
		// checks to see if the current line is sim data header 2 and if it 
		//is it will stop looping through each line
		if (strstr(data, SIMULATION_DATA_FILE_HEADER2) != NULL) break;
		// if the data array is not a NULL character, it will increment up after each line is read
		if (data != NULL) N++;
	}
	*NS = N; // pointer NS is set to line counter N
	printf("Number of cases = %d", *NS);
	rewind(fin); // rewinds the document from the beginning to read and tokenize
	// dynamically allocates memory for a sim data array 
	//to store simulation data and intializes it all to 0
	SD = (SIMULATION_DATA*)calloc(*NS, sizeof(SIMULATION_DATA));
	if (SD == NULL) exit(0); // if the array contains NULL, the program will end upon exit
	while (fgets(data, MAX_BUFF_SIZE, fin) != NULL) // reads through the input file line by line
	{
		// if the current line contains the sim data header string, it will skip
		if (strstr(data, SIMULATION_DATA_FILE_HEADER1)) continue;
		// if the current line contains dashes, it will break out of the loop
		if (strstr(data, DASHES)) break; 
	}
	// looping through for total number of simulations
	for (n = 0; n < *NS; n++)
	{
		// reads each line of the input file
		pRet = fgets(data, MAX_BUFF_SIZE, fin);
		// first token of the line
		tok = strtok_s(data, caseSeps, &nextToken); 
		// stores the first token of the case name into the sim data array strCase variable
		strcpy_s(SD[n].strCase, MAX_CASE_NAME_SIZE, tok);
		// tokenizes the next string and stores it into tok
		tok = strtok_s(NULL, caseSeps, &nextToken);
		// takes the second token being the width and 
		//converts the string into a double and stores into width variable
		SD[n].w = strtod(tok, &pGarbage); 
		// tokenizes the next string which is the height of the plate in the case
		tok = strtok_s(NULL, caseSeps, &nextToken);
		// converts the string into a double and stores 
		//it into the height variable of the sim data array
		SD[n].h = strtod(tok, &pGarbage); 
		// tokenizes the third string dx
		tok = strtok_s(NULL, caseSeps, &nextToken);
		// converts the string into a double and stores the 
		//dx value into the dx variable in sim data array
		SD[n].dx = strtod(tok, &pGarbage); 
		// tokenizes the final variable dy
		tok = strtok_s(NULL, caseSeps, &nextToken); 
		// converts the string into a double and stores the dy value 
		//into the dy variable in sim data array
		SD[n].dy = strtod(tok, &pGarbage);
	}
	// reads through each line of the input file
	while (fgets(data, MAX_BUFF_SIZE, fin) != NULL) 
	{
		// reading through each line of input file
		while (fgets(data, MAX_BUFF_SIZE, fin) != NULL) 
		{
			// checks the return of each line 
			pRet = fgets(data, MAX_BUFF_SIZE, fin);
			// if the current line is a blank line it'll skip past it
			if (isBlankLine(data)) continue;
			// if the current line is the header 2, it'll skip past
			if (strstr(data, SIMULATION_DATA_FILE_HEADER2) != NULL) continue;
			// if the current line is dashes it will break out of the loop
			if (strstr(data, DASHES) != NULL) break; 
		}
		// loops through for the total number of simulation cases
		for (i = 0; i < *NS; i++) 
		{
			// loop while the tokenized string is not the same as the 
			//ith element of the strcase and i is less than total number of sims
			while (tok != SD[i].strCase && i < *NS) 
			{
				// if the tokenized line is blank, it will break out of the loop
				if (isBlankLine(tok) != NULL) break;
				// reads the lines of the input file
				pRet = fgets(data, MAX_BUFF_SIZE, fin); 
				// intializes the n counter to 0
				n = 0; 
				// nested loop looping through the lines as long as they aren't
				//null and as long as n is less than the total number of walls
				while (fgets(data, MAX_BUFF_SIZE, fin) != NULL && n < NUM_WALLS) 
				{
					// tokenizes the first string in the case data line
					tok = strtok_s(data, caseSeps, &nextToken); 
					// tokenizes the next string in the cases data line
					tok = strtok_s(NULL, caseSeps, &nextToken); 
					// stores the case type as an integer that is created by the case to int function
					SD[i].bc[n].nType = caseTypetoInt(tok);
					// if the nType is equal to the corresponding int for type CONST, enter the loop
					if (SD[i].bc[n].nType == BC_TYPE_CONST) 
					{
						// tokenizes the subsequent string on the line 
						tok = strtok_s(NULL, caseSeps, &nextToken);
						// stores the tokenized string into Ta
						SD[i].bc[n].Ta = strtod(tok, &pGarbage); 
						tok = strtok_s(NULL, caseSeps, &nextToken);
						// stores the next token into za
						SD[i].bc[n].za = strtod(tok, &pGarbage);
						// tokenizes
						tok = strtok_s(NULL, caseSeps, &nextToken);
						// stores the next token into zb
						SD[i].bc[n].zb = strtod(tok, &pGarbage); 
					}
					// if the nType int is equal to int for type COSINE, enter the loop
					else if (SD[i].bc[n].nType == BC_TYPE_COSINE) 
					{
						// tokenizes
						tok = strtok_s(NULL, caseSeps, &nextToken); 
						// stores the token into Ta
						SD[i].bc[n].Ta = strtod(tok, &pGarbage);
						// tokenizes
						tok = strtok_s(NULL, caseSeps, &nextToken); 
						// stores teh token into za
						SD[i].bc[n].za = strtod(tok, &pGarbage); 
						// tokenizes
						tok = strtok_s(NULL, caseSeps, &nextToken);
						// stores the token into zb
						SD[i].bc[n].zb = strtod(tok, &pGarbage); 
					}
					// checks to see if the nType int is equal to int for type INSULATED
					else if (SD[i].bc[n].nType == BC_TYPE_INSULATED) 
					{
						// tokenizes
						tok = strtok_s(NULL, caseSeps, &nextToken); 
						// stores token into za
						SD[i].bc[n].za = strtod(tok, &pGarbage);
						// tokenizes
						tok = strtok_s(NULL, caseSeps, &nextToken); 
						// stores token into zb
						SD[i].bc[n].zb = strtod(tok, &pGarbage); 

					}
					// checks to see if the nType int is equal to int for type POLY
					else if (SD[i].bc[n].nType == BC_TYPE_POLY) 
					{
						// tokenizes
						tok = strtok_s(NULL, caseSeps, &nextToken);
						// stores token into Ta
						SD[i].bc[n].Ta = strtod(tok, &pGarbage); 
						// tokenizes
						tok = strtok_s(NULL, caseSeps, &nextToken); 
						// stores token into Tb
						SD[i].bc[n].Tb = strtod(tok, &pGarbage); 
						// tokenizes
						tok = strtok_s(NULL, caseSeps, &nextToken); 
						// stores token into za
						SD[i].bc[n].za = strtod(tok, &pGarbage); 
						// tokenizes
						tok = strtok_s(NULL, caseSeps, &nextToken);
						// stores token into zb
						SD[i].bc[n].zb = strtod(tok, &pGarbage);
						// tokenizes
						tok = strtok_s(NULL, caseSeps, &nextToken); 
						// stores token into ma
						SD[i].bc[n].ma = strtod(tok, &pGarbage);
						// tokenizes
						tok = strtok_s(NULL, caseSeps, &nextToken); 
						// stores token into mb
						SD[i].bc[n].mb = strtod(tok, &pGarbage); 
					}
					// checks to see if nType is equal to the int for type SINE
					else if (SD[i].bc[n].nType == BC_TYPE_SINE) 
					{
						// tokenizes
						tok = strtok_s(NULL, caseSeps, &nextToken);
						// stores token into Ta
						SD[i].bc[n].Ta = strtod(tok, &pGarbage);
						// tokenizes
						tok = strtok_s(NULL, caseSeps, &nextToken);
						// stores token into k
						SD[i].bc[n].k = strtod(tok, &pGarbage);
						// tokenizes
						tok = strtok_s(NULL, caseSeps, &nextToken);
						// stores token into za
						SD[i].bc[n].za = strtod(tok, &pGarbage);
						// tokenizes
						tok = strtok_s(NULL, caseSeps, &nextToken);
						// stores token into zb
						SD[i].bc[n].zb = strtod(tok, &pGarbage); 
					}
					else break; // breaks out of the loop after each wall has been determined
					n++; // counts up the wall count n

				}
				i++; // counts up the i until it becomes greater than NS
			}
			break; // breaks out of the loop
			rewind(fin); // rewinds for each loop
		}
	}
	return SD;
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Prompts the user to make a selection based on which case they want to run the sim for 
//               and stores the selection that they made for use in the SD array
// ARGUMENTS:    P:  the 2D PLATEPOINT array
//               SD: the simulation data for case B
// RETURN VALUE: iS: 
int getUserSimulationChoice(SIMULATION_DATA* SD, int NS)
{
	int n; // loop counter
	bool bHasGarbage = false; // initilized boolean variable to false to check for garbage
	int ret, iS; // return and sim index variables

	system("cls"); // clears screen for the menu to print
	printHorizontalBorder(TL, TR); // prints the outside top border for the menu
	printf("%c", VL); // prints left vertical line
	printf("   %s   ", SIMULATION_CHOICE); // prints the simulation choice header
	printf("%c\n", VL); // right vertical line
	for (n = 0; n < NS; n++) drawStringLine(SD[n].strCase, n); // prints each simulation case choice
	printHorizontalBorder(BL, BR); // bottom menu border
	printf("\nSelect Simulation Case: "); // prompt to make simulation choice based on options 
	ret = scanf_s("%d", &iS); // checks to see if the return is only 1
	if (ret != 1) EXIT_FAILURE; // if the user enters more than one value, it won't work
	if (iS != '\n') bHasGarbage = flushInputBuffer2(); // checks to see if the input is a garbage value
	// if the user choice is greater than the total number of sims avaiable or less than 0, prompt to try again
	if (iS > NS || iS <= 0) 
	{
		printf("Try making another choice!"); // Message to try again
		waitForEnterKey(); // holds the menu open until the ENTER key is pressed
	}
	return iS - 1; // returns the user menu selection minus one to account for the array starting at 0
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Dynamically allocates the memory for each array and intializes to 0 for each element 
//               and it intializes each x and y position of each node
// ARGUMENTS:    P:  the 2D PLATEPOINT array
//               SD: the simulation data for case B
// RETURN VALUE: PLATEPOINT P
PLATEPOINT** initialize(int iS, SIMULATION_DATA* SD, PLATEPOINT** P)
{
	double  w = SD[iS].w;//auxiliary variables for simulation data variables
	double dx = SD[iS].dx;
	double  h = SD[iS].h;
	double dy = SD[iS].dy;
	size_t n, i, j;
	//calculating the number of nodes in I and J
	SD[iS].I = nint((w / dx) + 1.0);
	SD[iS].J = nint((h / dy) + 1.0);
	//assigning I and J variables
	size_t I = SD[iS].I;
	size_t J = SD[iS].J;
	//1D PLATEPOINT array is allocated for size of I
	P = (PLATEPOINT**)calloc(I, sizeof(PLATEPOINT*));
	//If P is NULL exit the program
	if (P == NULL) exit(0);
	//allocate a 2D PLATEPOINT array
	for (i = 0; i < I; i++)
	{
		//allocate a 2D PLATEPOINT array for size of J
		P[i] = (PLATEPOINT*)calloc(J, sizeof(PLATEPOINT));
		//If P is NULL exit the program
		if (P == NULL) exit(0);
	}
	//If I or J are 0 exit the program
	if (I == 0 || J == 0) exit(0);
	//loop both i and j for the sizes I and J
	for (i = 0; i < I; i++)
	{
		for (j = 0; j < J; j++)
		{	//initialize x and y for the sizes dx and dy
			P[i][j].x = (double)i * dx;
			P[i][j].y = (double)j * dy;
		}
	}

	return P;
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Sets the boundary conditions for each case type depending on the users selection 
//               Sets Temps, za, zb, k, and other variables depending on which case type it is
// ARGUMENTS:    P:  the 2D PLATEPOINT array
//               iS: the user simulation selection
//               SD: the simulation data for the selected case
// RETURN VALUE: P
PLATEPOINT** SetBoundaryConditions(PLATEPOINT** P, SIMULATION_DATA* SD, int iS)
{
	size_t n, i = 0, j = 0; // counter variables
	size_t I = SD[iS].I, J = SD[iS].J; // auxilary variables for total number of nodes I and J
	// auxilary variables
	double w = SD[iS].w;
	double dx = SD[iS].dx;
	double h = SD[iS].h;
	double dy = SD[iS].dy;
	double T;

	// looping through to set the boundary conditions along each wall
	for (n = 0; n < NUM_WALLS; n++)
	{
		// check to see if the wall is type constant
		if (SD[iS].bc[n].nType == BC_TYPE_CONST)
		{
			if (n == TOP) // if the wall is the top wall
			{
				// auxilary variables
				double Tc = SD[iS].bc[n].Ta;
				double xa = SD[iS].bc[n].za;
				double xb = SD[iS].bc[n].zb;
				for (i = 0; i < I; i++) // looping through each node to set the temperatures
				{
					// if the current x is less than the xa or greater than the xb 
					//values based on the case, it will set the temp of the node = 0
					if (P[i][J - 1].x <= xa || P[i][J - 1].x >= xb)
					{
						T = 0;
						// setting both T_fd and T_a to equal the temp T
						P[i][J - 1].T_a = T;
						P[i][J - 1].T_fd = T;
					}
					if (P[i][J - 1].x >= xa && P[i][J - 1].x <= xb)
					{
						// sets both T_a and T_fd to the temperature of that case
						P[i][J - 1].T_a = Tc;
						P[i][J - 1].T_fd = P[i][J - 1].T_a;
					}
				}
			}
			else if (n == BOTTOM) // if the current wall is the bottom wall
			{
				// auxilary variables
				double Tc = SD[iS].bc[n].Ta;
				double xa = SD[iS].bc[n].za;
				double xb = SD[iS].bc[n].zb;
				double T;
				// looping through each node to set the conditions
				for (i = 0; i < I; i++)
				{
					double T = P[i][0].T_fd;

					// if the current x is less than the xa or greater than the xb 
					//values based on the case, it will set the temp of the node = 0
					if (P[i][0].x <= xa || P[i][0].x >= xb)
					{
						// setting both T_fd and T_a to equal the temp T
						T = 0;
						P[i][0].T_a = T;
						P[i][0].T_fd = T;
					}
					if (P[i][0].x >= xa && P[i][0].x <= xb)
					{
						// sets both T_a and T_fd to the temperature of that case
						P[i][0].T_a = Tc;
						P[i][0].T_fd = P[i][0].T_a;
					}
				}
			}
			else if (n == LEFT) // if the current wall is the left wall
			{
				// auxilary variables
				double Tc = SD[iS].bc[n].Ta;
				double ya = SD[iS].bc[n].za;
				double yb = SD[iS].bc[n].zb;
				double T;
				// looping through the J nodes for the j direction to set the temps
				for (j = 0; j < J; j++)
				{
					double T = P[0][j].T_fd;
					// if the current y is less than the ya or greater than the yb 
					//values based on the case, it will set the temp of the node = 0
					if (P[0][j].y <= ya || P[0][j].y >= yb)
					{
						// setting both T_fd and T_a to equal the temp T
						T = 0;
						P[0][j].T_a = T;
						P[0][j].T_fd = T;
					}
					if (P[0][j].y >= ya && P[0][j].y <= yb)
					{
						// sets both T_a and T_fd to the temperature of that case
						P[0][j].T_a = Tc;
						P[0][j].T_fd = P[0][j].T_a;
					}
				}
			}
			else if (n == RIGHT) // if the current wall is the right wall
			{
				// auxilary variable
				double Tc = SD[iS].bc[n].Ta;
				double ya = SD[iS].bc[n].za;
				double yb = SD[iS].bc[n].zb;
				double T;
				// looping through the J nodes for the j direction to set the temps
				for (j = 0; j < J; j++)
				{
					double T = P[I - 1][j].T_fd;

					// if the current y is less than the ya or greater than the yb 
					//values based on the case, it will set the temp of the node = 0
					if (P[I - 1][j].y <= ya || P[0][j].y >= yb)
					{
						// setting both T_fd and T_a to equal the temp T
						T = 0;
						P[I - 1][j].T_a = T;
						P[I - 1][j].T_fd = T;
					}
					if (P[I - 1][j].y >= ya && P[0][j].y <= yb)
					{
						// sets both T_a and T_fd to the temperature of that case
						P[I - 1][j].T_a = Tc;
						P[I - 1][j].T_fd = Tc;
					}
				}
			}
		}
		else if (SD[iS].bc[n].nType == BC_TYPE_COSINE) // if the case type is cosine for the wall
		{
			if (n == TOP) // if it is a top wall
			{
				// auxilary variables
				double Tm = SD[iS].bc[n].Ta;
				double xa = SD[iS].bc[n].za;
				double xb = SD[iS].bc[n].zb;
				double T;
				
				// loops through the i nodes to set the temps
				for (i = 0; i < I; i++)
				{
					// if the current x is less than the xa or greater than the xb 
					//values based on the case, it will set the temp of the node = 0
					if (P[i][J - 1].x <= xa || P[i][J - 1].x >= xb)
					{
						// setting both T_fd and T_a to equal the temp T
						T = 0;
						P[i][J - 1].T_a = T;
						P[i][J - 1].T_fd = T;
					}
					if(P[i][J - 1].x >= xa && P[i][J - 1].x <= xb)
					{
						// uses the cosine function to calculate the temps 
						//and sets the nodes to that temp for T_fd and T_a
						T = (Tm / 2.0) * (1.0 - cos(2.0 * PI * ((P[i][J - 1].x - xa) / (xb - xa))));
						P[i][J - 1].T_a = T;
						P[i][J - 1].T_fd = T;
					}
				}
			}
			else if (n == BOTTOM)
			{
				// auxilary variables
				double Tm = SD[iS].bc[n].Ta;
				double xa = SD[iS].bc[n].za;
				double xb = SD[iS].bc[n].zb;
				double T;

				// loops through the i nodes to set the temps
				for (i = 0; i < I; i++)
				{
					// if the current x is less than the xa or greater than the xb 
					//values based on the case, it will set the temp of the node = 0
					if (P[i][0].x <= xa || P[i][0].x >= xb)
					{
						// setting both T_fd and T_a to equal the temp T
						T = 0;
						P[i][0].T_a = T;
						P[i][0].T_fd = T;
					}
					if (P[i][0].x >= xa && P[i][0].x <= xb)
					{
						// uses the cosine function to calculate the temps and 
						//sets the nodes to that temp for T_fd and T_a
						T = (Tm / 2.0) * (1.0 - cos(2.0 * PI * ((P[i][0].x - xa) / (xb - xa))));
						P[i][0].T_a = T;
						P[i][0].T_fd = T;
					}
				}
			}
			else if (n == LEFT) // if the current wall is the left wall
			{
				// auxilary variables
				double Tm = SD[iS].bc[n].Ta;
				double ya = SD[iS].bc[n].za;
				double yb = SD[iS].bc[n].zb;
				double T;

				// loops through the j nodes to set the temps
				for (j = 0; j < J; j++)
				{
					// if the current y is less than the ya or greater than the yb 
					//values based on the case, it will set the temp of the node = 0
					if (P[0][j].y <= ya || P[0][j].y >= yb)
					{
						// setting both T_fd and T_a to equal the temp T
						T = 0;
						P[0][j].T_a = T;
						P[0][j].T_fd = T;
					}
					if (P[0][j].y >= ya && P[0][j].y <= yb)
					{
						// uses the cosine function to calculate the temps and sets 
						//the nodes to that temp for T_fd and T_a
						T = (Tm / 2.0) * (1.0 - cos(2.0 * PI * ((P[0][j].y - ya) / (yb - ya))));
						P[0][j].T_a = T;
						P[0][j].T_fd = P[0][j].T_a;
					}
				}
			}
			else if (n == RIGHT)
			{
				// auxilary variables
				double Tm = SD[iS].bc[n].Ta;
				double ya = SD[iS].bc[n].za;
				double yb = SD[iS].bc[n].zb;
				double T;

				// loops through the j nodes to set the temps
				for (j = 0; j < J; j++)
				{
					// if the current y is less than the ya or greater than the yb values 
					//based on the case, it will set the temp of the node = 0
					if (P[I - 1][j].y <= ya || P[I - 1][j].y >= yb)
					{
						// setting both T_fd and T_a to equal the temp T
						T = 0;
						P[I - 1][j].T_a = T;
						P[I - 1][j].T_fd = T;
					}
					if (P[I - 1][j].y >= ya && P[I - 1][j].y <= yb)
					{
						// uses the cosine function to calculate the temps and sets the nodes 
						//to that temp for T_fd and T_a
						T = (Tm / 2.0) * (1.0 - cos(2.0 * PI * ((P[I - 1][j].y - ya) / (yb - ya))));
						P[I - 1][j].T_a = T;
						P[I - 1][j].T_fd = T;
					}
				}
			}
		}
		else if (SD[iS].bc[n].nType == BC_TYPE_INSULATED) // if the wall type is type insulated
		{
			double T;

			// loops through the j nodes to set the temps
			for (j = 1; j < J - 1; j++)
			{
				// setting both T_fd and T_a to equal the temp T
				T = 0;
				P[I - 1][j].T_a = T;
				P[I - 1][j].T_fd = T;
			}
		}
		else if (SD[iS].bc[n].nType == BC_TYPE_POLY) // if the wall type is type poly
		{
			// Auxilary variables
			double Ta = SD[iS].bc[n].Ta;
			double Tb = SD[iS].bc[n].Tb;
			double ma = SD[iS].bc[n].ma;
			double mb = SD[iS].bc[n].mb;

			if (n == TOP) // if the wall is the top wall
			{
				// auxilary variable
				double  xa = SD[iS].bc[n].za;
				double  xb = SD[iS].bc[n].zb;
				double del = xb - xa;
				double   a = Ta;
				double   b = ma * del;
				double   c = 3.0 * (Tb - Ta) - (2.0 * ma + mb) * del;
				double   d = -2.0 * (Tb - Ta) + (ma + mb) * del;
				double T;

				// loops through the i nodes to set the temps
				for (i = 0; i < I; i++)
				{
					// if the current phi is less than the ya or greater than the yb values 
					//based on the case, it will set the temp of the node = 0
					double phi = (P[i][J - 1].x - xa) / del;
					if (phi >= 0.0 && phi <= 1.0)
					{
						// calculates the temp of each node using a polynomial equation
						T = a + b * phi + c * phi * phi + d * phi * phi * phi;
						P[i][J - 1].T_a = T;
						P[i][J - 1].T_fd = T;
					}
					if (phi <= 0.0 || phi >= 1.0)
					{
						// setting both T_fd and T_a to equal the temp T
						T = 0;
						P[i][J - 1].T_a = T;
						P[i][J - 1].T_fd = T;
					}
				}
			}
			else if (n == BOTTOM) // if the current wall is the bottom wall
			{
				// auxilary variable
				double  xa = SD[iS].bc[n].za;
				double  xb = SD[iS].bc[n].zb;
				double del = xb - xa;
				double   a = Ta;
				double   b = ma * del;
				double   c = 3.0 * (Tb - Ta) - (2.0 * ma + mb) * del;
				double   d = -2.0 * (Tb - Ta) + (ma + mb) * del;
				double T;

				// loops through the i nodes to set the temps
				for (i = 0; i < I; i++)
				{
					// if the current phi is less than the ya or greater than the yb 
					//values based on the case, it will set the temp of the node = 0
					double phi = (P[i][0].x - xa) / del;
					if (phi >= 0.0 && phi <= 1.0)
					{
						// calculates the temp of each node using a polynomial equation
						T = a + b * phi + c * phi * phi + d * phi * phi * phi;
						P[i][0].T_a = T;
						P[i][0].T_fd = T;
					}
					if (phi <= 0.0 || phi >= 1.0)
					{
						// setting both T_fd and T_a to equal the temp T
						T = 0;
						P[i][0].T_a = T;
						P[i][0].T_fd = T;
					}
				}
			}
			else if (n == LEFT) // if the current wall is the left wall
			{
				// auxilar variables
				double  ya = SD[iS].bc[n].za;
				double  yb = SD[iS].bc[n].zb;
				double del = yb - ya;
				double   a = Ta;
				double   b = ma * del;
				double   c = 3.0 * (Tb - Ta) - (2.0 * ma + mb) * del;
				double   d = -2.0 * (Tb - Ta) + (ma + mb) * del;
				double T;

				// loops through the j nodes to set the temps
				for (j = 0; j < J; j++)
				{
					// if the current phi is less than the ya or greater than the yb 
					//values based on the case, it will set the temp of the node = 0
					double phi = (P[i][0].y - ya) / del;
					if (phi >= 0.0 && phi <= 1.0)
					{
						// calculates the temp of each node using a polynomial equation
						T = a + b * phi + c * phi * phi + d * phi * phi * phi;
						P[0][j].T_a = T;
						P[0][j].T_fd = T;
					}
					if (phi <= 0.0 || phi >= 1.0)
					{
						// setting both T_fd and T_a to equal the temp T
						T = 0;
						P[0][j].T_a = T;
						P[0][j].T_fd = T;
					}
				}
			}
			else if (n == RIGHT) // if the current wall is the right wall
			{
			    // auxilary variables
				double  ya = SD[iS].bc[n].za;
				double  yb = SD[iS].bc[n].zb;
				double del = yb - ya;
				double   a = Ta;
				double   b = ma * del;
				double   c = 3.0 * (Tb - Ta) - (2.0 * ma + mb) * del;
				double   d = -2.0 * (Tb - Ta) + (ma + mb) * del;
				double T;

				// loops through the j nodes to set the temps
				for (j = 0; j < J; j++)
				{
					// if the current phi is less than the ya or greater than the yb 
					//values based on the case, it will set the temp of the node = 0
					double phi = (P[i][J - 1].y - ya) / del;
					if (phi >= 0.0 && phi <= 1.0)
					{
						// calculates the temp of each node using a polynomial equation
						T = a + b * phi + c * phi * phi + d * phi * phi * phi;
						P[I - 1][j].T_a = T;
						P[I - 1][j].T_fd = T;
					}
					if (phi <= 0.0 || phi >= 1.0)
					{
						// setting both T_fd and T_a to equal the temp T
						T = 0;
						P[I - 1][j].T_a = T;
						P[I - 1][j].T_fd = T;
					}
				}
			}
		}
		else if (SD[iS].bc[n].nType == BC_TYPE_SINE) // if the wall is type sine
		{
			if (n == TOP) // if the current wall is a top wall
			{
				// auxilary variables
				double Ta = SD[iS].bc[n].Ta;
				double xa = SD[iS].bc[n].za;
				double xb = SD[iS].bc[n].zb;
				double k = SD[iS].bc[n].k;
				double T;
				// looping through the i nodes
				for (i = 0; i < I; i++)
				{
					// checks to see if the current x value is greater than xb or less than xa
					if (P[i][J - 1].x <= xa || P[i][J - 1].x >= xb)
					{
						// setting both T_fd and T_a to equal the temp T
						T = 0;
						P[i][J - 1].T_a = T;
						P[i][J - 1].T_fd = T;
					}
					if (P[i][J - 1].x >= xa && P[i][J - 1].x <= xb)
					{
						// uses the sinusoidal function to calculate the temps and 
						//sets the nodes to that temp for T_fd and T_a
						T = Ta * sin(k * PI * ((P[i][j].x - xa) / (xb - xa)));
						P[i][J - 1].T_a = T;
						P[i][J - 1].T_fd = T;
					}
				}
			}
			else if (n == BOTTOM)
			{
				double Ta = SD[iS].bc[n].Ta;
				double xa = SD[iS].bc[n].za;
				double xb = SD[iS].bc[n].zb;
				double k = SD[iS].bc[n].k;
				double T;

				// looping through the i nodes
				for (i = 0; i < I; i++)
				{
					// checks to see if the current x value is greater than xb or less than xa
					if (P[i][0].x <= xa || P[i][0].x >= xb)
					{
						// setting both T_fd and T_a to equal the temp T
						T = 0;
						P[i][0].T_a = T;
						P[i][0].T_fd = T;
					}
					if (P[i][0].x >= xa && P[i][0].x <= xb)
					{
						// uses the sinusoidal function to calculate the temps and sets the 
						//nodes to that temp for T_fd and T_a
						T = Ta * sin(k * PI * ((P[i][0].x - xa) / (xb - xa)));
						P[i][0].T_a = T;
						P[i][0].T_fd = T;
					}
				}
			}
			else if (n == LEFT) // if the current wall is the left wall
			{
				// auxilary variables
				double Ta = SD[iS].bc[n].Ta;
				double ya = SD[iS].bc[n].za;
				double yb = SD[iS].bc[n].zb;
				double k = SD[iS].bc[n].k;
				double T;

				// looping through the j nodes
				for (j = 0; j < J; j++)
				{
					// checks to see if the current y value is greater than yb or less than ya
					if (P[0][j].y <= ya || P[0][j].y >= yb)
					{
						// setting both T_fd and T_a to equal the temp T
						T = 0;
						P[0][j].T_a = T;
						P[0][j].T_fd = T;
					}
					if (P[0][j].y >= ya && P[0][j].y <= yb)
					{
						// uses the sinusoidal function to calculate the temps 
						//and sets the nodes to that temp for T_fd and T_a
						T = Ta * sin(k * PI * ((P[0][j].y - ya) / (yb - ya)));
						P[0][j].T_a = T;
						P[0][J].T_fd = T;
					}
				}
			}
			else if (n == RIGHT) // if the current wall is the right wall
			{
				// auxilary variables
				double Ta = SD[iS].bc[n].Ta;
				double ya = SD[iS].bc[n].za;
				double yb = SD[iS].bc[n].zb;
				double k = SD[iS].bc[n].k;
				double T;

				// looping through the j nodes
				for (j = 0; j < J; j++)
				{
					// checks to see if the current y value is greater than yb or less than ya
					if (P[I - 1][j].y <= ya || P[I - 1][j].y >= yb)
					{
						// setting both T_fd and T_a to equal the temp T
						T = 0;
						P[I - 1][j].T_a = T;
						P[I - 1][j].T_fd = T;
					}
					if (P[I - 1][j].y >= ya && P[I - 1][j].y <= yb)
					{
						// uses the sinusoidal function to calculate the temps and sets 
						//the nodes to that temp for T_fd and T_a
						T = Ta * sin(k * PI * ((P[I - 1][j].y - ya) / (yb - ya)));
						P[I - 1][j].T_a = T;
						P[I - 1][j].T_fd = T;
					}
				}
			}
		}
	}
	// calculates the average temparture of the top left node
	P[0][J - 1].T_a = (P[0][J - 2].T_a + P[1][J - 1].T_a) / 2;
	P[0][J - 1].T_fd = P[0][J - 1].T_a;

	// calculates the average temperature of the bottom left node
	P[0][0].T_a = (P[0][1].T_a + P[1][0].T_a) / 2.0;
	P[0][0].T_fd = P[0][0].T_a;
	// if the case is insulated, it won't take the average temperature of the top right node
	if (SD[iS].bc[RIGHT].nType != BC_TYPE_INSULATED)
	{
		// top right node
		P[I - 1][J - 1].T_a = (P[I - 2][J - 1].T_a + P[I - 1][J - 2].T_a) / 2.0;
		P[I - 1][J - 1].T_fd = P[I - 1][J - 1].T_a;

		// bottom right node
		P[I - 1][0].T_a = (P[I - 2][0].T_a + P[I - 1][1].T_a) / 2.0;
		P[I - 1][0].T_fd = P[I - 1][0].T_a;
	}

	return P;
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Uses Finite-difference method to numerically solve for the temperature of each node 
//               Cycles through each node and finds the temperature based on the average of neighbouring nodes
// ARGUMENTS:    P:  the 2D PLATEPOINT array
//               SD: the simulation data for the selected case
// RETURN VALUE: none
void GetNumericalSolution(PLATEPOINT** P, const SIMULATION_DATA SD)
{
	FILE* fConverge = NULL;
	errno_t err;
	double rmax = 0; // defines and initializes rmax to zero
	int I = nint((SD.w / SD.dx) + 1.0); // calculates values for I
	int J = nint((SD.h / SD.dy) + 1.0); // calculates values for J
	int i = 0, j = 0; // counters 
	char strConvergenceFile[MAX_BUFF_SIZE]; // convergence file string name
	double RMS = 0.0; // variable holder for RMS value
	double dy = SD.dy; // defines dy from structure 
	double dx = SD.dx; // defines dx from structure 
	double lamda = pow(SD.dx / SD.dy, 2.0); // calculates lamda 
	int iter = 0; // iteration counter

	sprintf_s(strConvergenceFile, MAX_BUFF_SIZE, "%s convergence.dat", SD.strCase);
	err = fopen_s(&fConverge, strConvergenceFile, "w");
	if (err != 0 || fConverge == NULL)
	{
		printf("Cannot open \"%s\" for writing...", strConvergenceFile);
		waitForEnterKey();
		EXIT_FAILURE;
	}

	do
	{
		for (j = 1; j < J - 1; j++) // sweeping through the nodes vertically 
		{
			for (i = 1; i < I - 1; i++) // sweeping through the nodes horizontally 
			{
				// calculates value for temperature finite difference by using the formula found in 
				//finite difference laplace.pdf
				P[i][j].T_fd = (P[i + 1][j].T_fd + P[i - 1][j].T_fd + lamda * (P[i][j + 1].T_fd + 
					P[i][j - 1].T_fd)) / (2.0 * (1.0 + lamda));
			}
			if (SD.bc[RIGHT].nType == BC_TYPE_INSULATED) // special formula used for insulated right wall
			{
				P[I - 1][j].T_fd = (2.0 * P[I - 2][j].T_fd + lamda * (P[I - 1][j + 1].T_fd + 
					P[I - 1][j - 1].T_fd)) / (2.0 * (1.0 + lamda));
			}
		}
		RMS = 0.0; // resets RMS to zero
		rmax = 0.0; // resets rmax to zero
		for (j = 1; j < J - 1; j++) // sweeping through the nodes vertically 
		{
			for (i = 1; i < I - 1; i++) // sweeping through the nodes horizontally 
			{
				// calculates value for residual by using the formula found in finite difference laplace.pdf
				P[i][j].res = fabs(P[i][j].T_fd - (P[i + 1][j].T_fd + P[i - 1][j].T_fd + lamda * 
					(P[i][j + 1].T_fd + P[i][j - 1].T_fd)) / (2.0 * (1.0 + lamda)));
				// if the resiudal is greater than the current rmax, replace the rmax with residual
				if (P[i][j].res > rmax) rmax = P[i][j].res;
				// add the calculated value onto the previous value 
				RMS += pow(P[i][j].res, 2.0); 
			}
			if (SD.bc[RIGHT].nType == BC_TYPE_INSULATED) // for insulated right wall 
			{
				P[I - 1][j].res = fabs(P[i][j].T_fd - (2.0 * P[I - 2][j].T_fd + lamda * 
					(P[I - 1][j + 1].T_fd + P[I - 1][j - 1].T_fd)) / (2.0 * (1.0 + lamda)));
				if (P[i][j].res > rmax) rmax = P[i][j].res;
				RMS += pow(P[i][j].res, 2.0);
			}
		}
		RMS = sqrt(RMS / (((double)I - 2) * ((double)J - 2))); // calculates RMS
		iter++; // iter increments 
	    // do the loop while iter is less than or equal to MAX_ITER AND rmax is 
		//greater or eqal to MAX_RESIDUAL AND RMS greater or equal to MAX_RESIDUAL  
		if (fConverge == 0) exit(0);
		fprintf(fConverge, "%12.5le, %12.5le, %d\n", rmax, RMS, iter);

	} while (iter <= MAX_ITER && (rmax >= MAX_RESIDUAL && RMS >= MAX_RESIDUAL));

	// prints to screen - the values of iter, rmax and RMS
	printf("\nNumber of iterations: %d", iter);
	printf("\nRmax = %.5le", rmax);
	printf("\nRMS = %.5le\n\n", RMS);

	fclose(fConverge);
	printf("\nPrinted data to file \"%s\n", strConvergenceFile);
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Computes the analytical solution for case A and store the temperature values into the 
//               PLATEPOINT array.
// ARGUMENTS:    P:  the 2D PLATEPOINT array
//               SD: the simulation data for case B
// RETURN VALUE: none
void GetCaseAAnalyticalSolution(PLATEPOINT** P, const SIMULATION_DATA* SD)
{
	size_t i, j, n;                                       // loop counters
	size_t I = SD->I, J = SD->J;                       // 2D array dimensions for the simulation case
	double T, T1 = SD->bc[TOP].Ta;                     // temperature, peak temperature
	double h = SD->h, w = SD->w, k = SD->bc[TOP].k;    // plate height/width, k factor in sine function
	double x, y;                                       // plate coordinates for 2D array element

	for (i = 1; i < I - 1; i++) // boundaries already done!
	{
		for (j = 1; j < J - 1; j++)
		{
			// auxilary variables so Temperature calculation formula can be written on one line
			x = P[i][j].x;
			y = P[i][j].y;
			//reset Tsum to 0
			double Tsum = 0;
			//Iterate the infinite sum for 100 times
			for (n = 1; n < 100; n++)
			{
				double A;//used as auxiliary variables 
				double B;
				double C;
				double D;
				//Infinite sum formula is broken into smaller variables to add loop break condition for sinh
				A = (1 - cos(n * PI)) / n;
				B = sin(n * PI * x / w);
				C = sinh(n * PI * y / w);
				D = sinh(n * PI * h / w);
				//if sinh gets bigger than DBL_MAX, then break loop
				if (D > DBL_MAX) break;
				//infinite sum formula combined
				T = A * B * C / D;
				//sum the temperature value
				Tsum += T;
			}
			//Temperature formula for case 2: constant temperature on the upper body
			//After iterating the infinite sum, calculate the temperature
			T = T0 + 2 / PI * (T1 - T0) * Tsum;
			//Store temperature value to PLATEPOINT array
			P[i][j].T_a = T;
		}
	}
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Computes the analytical solution for case B and store the temperature values into the 
//               PLATEPOINT array.
// ARGUMENTS:    P:  the 2D PLATEPOINT array
//               SD: the simulation data for case B
// RETURN VALUE: none
void GetCaseBAnalyticalSolution(PLATEPOINT** P, const SIMULATION_DATA* SD)
{
	size_t i, j;                                       // loop counters
	size_t I = SD->I, J = SD->J;                     // 2D array dimensions for the simulation case
	double T, T1 = SD->bc[TOP].Ta;                    // temperature, peak temperature
	double h = SD->h, w = SD->w, k = SD->bc[TOP].k; // plate height/width, k factor in sine function
	double x, y;                                       // plate coordinates for 2D array element

	for (i = 1; i < I - 1; i++) // boundaries already done!
	{
		for (j = 1; j < J - 1; j++)
		{
			// auxilary variables so Temperature calculation formula can be written on one line
			x = P[i][j].x;
			y = P[i][j].y;
			//Temperature for case 1: sinusoidal distribution on the upper boundary
			T = T0 + T1 * sin(k * PI * x / w) * sinh(k * PI * y / w) / sinh(k * PI * h / w); // sine formula
			P[i][j].T_a = T;
		}
	}
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Computes the analytical solution for case C and store the temperature values into the 
//               PLATEPOINT array.
// ARGUMENTS:    P:  the 2D PLATEPOINT array
//               SD: the simulation data for case B
// RETURN VALUE: none
void GetCaseCAnalyticalSolution(PLATEPOINT** P, const SIMULATION_DATA* SD)
{
	size_t i, j;                                       // loop counters
	size_t I = SD->I, J = SD->J;                     // 2D array dimensions for the simulation case
	double T, T1 = SD->bc[TOP].Ta;                    // temperature, peak temperature
	double h = SD->h, w = SD->w, k = SD->bc[TOP].k; // plate height/width, k factor in sine function
	double x, y;                                       // plate coordinates for 2D array element

	for (i = 1; i < I; i++) // boundaries already done!
	{
		for (j = 1; j < J - 1; j++)
		{
			// auxilary variables so Temperature calculation formula can be written on one line
			x = P[i][j].x;
			y = P[i][j].y;
			//Temperature formula for case 3
			T = T0 + T1 * ((sin((k - 1 / 2) * PI * x / w) * sinh((k - 1 / 2) 
				* PI * y / w)) / sinh((k - 1 / 2) * PI * h / w));
			P[i][j].T_a = T;
		}
	}
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Prints the solution to Matlab in order to display and graph the temperature distribution 
//               onto the steel plate
// ARGUMENTS:    P:  the 2D PLATEPOINT array
//               SD: the simulation data for case B
// RETURN VALUE: none
void printSolution(PLATEPOINT** P, const SIMULATION_DATA* pSD)
{

	size_t i, j;                                 // loop/temp variables
	FILE* fa = NULL, * ffd = NULL, * fres = NULL;  // for printing solutions to file for matlab
	char strFileNameAnalytical[MAX_BUFF_SIZE];   // buffer to hold analytical output file name
	char strFileNameFD[MAX_BUFF_SIZE];           // buffer to hold F-D output file name
	char strFileNameResidual[MAX_BUFF_SIZE];     // buffer to hold Residual output file name
	errno_t err;                                 // check fopen

	// open analytical output file only if applicable!  create name dynamically with sprintf_s
	if (pSD->nCaseType != CASE_TYPE_TEST)
	{
		sprintf_s(strFileNameAnalytical, MAX_BUFF_SIZE, "%s Analytical.dat", pSD->strCase);
		err = fopen_s(&fa, strFileNameAnalytical, "w");
		if (err != 0 || fa == NULL)
		{
			printf("Cannot open \"%s\" for writing. Skipping printout...\n", strFileNameAnalytical);
			return;
		}
	}

	// open finite different output file.   create name dynamically with sprintf_s
	sprintf_s(strFileNameFD, MAX_BUFF_SIZE, "%s Finite Difference.dat", pSD->strCase);
	err = fopen_s(&ffd, strFileNameFD, "w");
	if (err != 0 || ffd == NULL)
	{
		printf("Cannot open \"%s\" for writing. Skipping printout...\n", strFileNameFD);
		return;
	}

	// open residual output file.  create name dynamically with sprintf_s
	sprintf_s(strFileNameResidual, MAX_BUFF_SIZE, "%s Residual.dat", pSD->strCase);
	err = fopen_s(&fres, strFileNameResidual, "w");
	if (err != 0 || fres == NULL)
	{
		printf("Cannot open \"%s\" for writing. Skipping printout...\n", strFileNameResidual);
		return;
	}

	//--------------  print outputs to dat files for matlab -------------------------
	for (i = 0; i < pSD->I; i++)
	{
		for (j = 0; j < pSD->J; j++)
		{
			fprintf(fres, "%+12.5le,%+12.5le,%+12.5le\n", P[i][j].x, P[i][j].y, P[i][j].res);
			fprintf(ffd, "%+12.5le,%+12.5le,%+12.5le\n", P[i][j].x, P[i][j].y, P[i][j].T_fd);

			// don't print analytical if there isn't a solution!
			if (pSD->nCaseType != CASE_TYPE_TEST)
				fprintf(fa, "%+12.5le,%+12.5le,%+12.5le\n", P[i][j].x, P[i][j].y, P[i][j].T_a);
		}
	}

	// close the files
	if (pSD->nCaseType != CASE_TYPE_TEST) fclose(fa);
	fclose(ffd);
	fclose(fres);

	// echo success to screen
	if (pSD->nCaseType != CASE_TYPE_TEST) printf("Printed data to \"%s\"\n", strFileNameAnalytical);
	printf("Printed data to \"%s\"\n", strFileNameFD);
	printf("Printed data to \"%s\"\n", strFileNameResidual);
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Frees the memory that was allocated to the dynamic arrays for the platepoint and SD struc 
// ARGUMENTS:    P:  the 2D PLATEPOINT array
//               SD: the simulation data for case B
//                I: total nodes in the i direction
// RETURN VALUE: none
void FreeMemory(PLATEPOINT** P, SIMULATION_DATA* SD, size_t I)
{
	size_t i;
	// Freeing SD Array
	free(SD);

	// Freeing Platepoint Array
	for (i = 0; i < I; i++)
		free(P[i]);
	free(P);

}


//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Computes the integer nearest to the given double.
// ARGUMENTS:    d:  the double
// RETURN VALUE: the nearest integer
int nint(double d)
{
	return (int)floor(d + 0.5);
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Checks if the string returned by fgets is a blank line (i.e., only whitespace characters)
// ARGUMENTS:    strLine:  the string
// RETURN VALUE: true if only whitespace found in string, false if not
bool isBlankLine(const char* strLine)
{
	char* tok, * nextTok = NULL;    // tokenize the string with whitespace to see if that's all there is
	const char* seps = " \t\n\r";
	char buff[MAX_BUFF_SIZE];      // copy the string because it's a const and strtok_s can't handle const

	strcpy_s(buff, MAX_BUFF_SIZE, strLine);
	tok = strtok_s(buff, seps, &nextTok);

	if (tok == NULL)
		return true;
	else
		return false;
}

//--------------------------------------------------------------------------------------------
// DESCRIPTION:  Prints a horizontal border
// ARGUMENTS:    cLeft: the left border character, cRight: rght border character
// RETURN VALUE: int numChars
int printHorizontalBorder(char cLeft, char cRight)  // prints a character repeatedly
{
	int n;  // loop counter
	int numChars = 0; // number of characters printed by printfs

	numChars += printf("%c", cLeft);  // left corner
	for (n = 0; n < strlen(SIMULATION_CHOICE) + 2 * TABLE_MARGIN_SIZE; n++) numChars += printf("%c", HL);  // mid section line
	numChars += printf("%c\n", cRight) - 1;  // right corner (-1 for newline)

	return numChars;
}

//--------------------------------------------------------------------------------------------
// DESCRIPTION:  Prints out a string 
// ARGUMENTS:    const char: string message
// RETURN VALUE: none
void drawStringLine(const char* string1, int n)
{
	int centright;
	int numString = strlen(string1);

	centright = (TABLE_COLUMN_WIDTH - TABLE_MARGIN_SIZE - numString);
	printf("%c", VL);
	printf("%c%c%c[%d]", ' ', ' ', ' ', n + 1);
	printf(" %-*s", strlen(SIMULATION_CHOICE) + TABLE_MARGIN_SIZE - 4, string1);
	printf("%c\n", VL);
}

//--------------------------------------------------------------------------------------------
// DESCRIPTION:  Converts the caseType string into an int value to compare with the 
//               nType integers
// ARGUMENTS:    char*: string message
// RETURN VALUE: int BCTYPE
int caseTypetoInt(char* string)
{
	int BCTYPE = 1;
	if (strcmp(string, "CONST") == 0)  BCTYPE = BC_TYPE_CONST;
	else if (strcmp(string, "COSINE") == 0) BCTYPE = BC_TYPE_COSINE;
	else if (strcmp(string, "INSULATED") == 0) BCTYPE = BC_TYPE_INSULATED;
	else if (strcmp(string, "POLY") == 0) BCTYPE = BC_TYPE_POLY;
	else if (strcmp(string, "SINE") == 0) BCTYPE = BC_TYPE_SINE;

	return BCTYPE;
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  This function flushes the input buffer to avoid scanf issues
//               ***** CALL THIS FUNCTION AFTER EVERY CALL TO SCANF!!! *****
// ARGUMENTS:    none
// RETURN VALUE: false if nothing or only '\n' is in the input buffer
//               true if extra keystrokes precede the '\n'.  Good for detecting left 
//               over garbage from scanf_s in the input buffer
bool flushInputBuffer2()
{
	unsigned char ch; // temp character variable
	bool bHasGarbage = false;

	// exit loop when all characters are flushed
	while ((ch = (unsigned char)getchar()) != '\n' && ch != EOF)
	{
		if (!bHasGarbage && !isspace(ch)) bHasGarbage = true;
	}
	return bHasGarbage;
}
//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Waits for user to press enter.  flushes stdin if keystroke precede enter
// ARGUMENTS:    none
// RETURN VALUE: none
void waitForEnterKey()
{
	unsigned char ch;
	if ((ch = (unsigned char)getchar()) != EOF && ch != '\n') flushInputBuffer2();
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Display a message and terminate the program
// ARGUMENTS:    the message string
// RETURN VALUE: none
void endProgram(const char* message)
{
	printf("\n");
	if (message != NULL) printf(message);
	printf("\nPress ENTER to end this program...");
	waitForEnterKey();
	exit(0);
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  print a character repeatedly
// ARGUMENTS:    ch: the character, N: the number of repititions
// RETURN VALUE: none
void printRepeatedChar(unsigned char ch, int N)
{
	int n;
	for (n = 0; n < N; n++) printf("%c", ch);
}

//------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Replaces the last character in a string, if it is a '\n', with '\0'.
// ARGUMENTS:    str:  the string
// RETURN VALUE: none
void removeNewline(char* str)
{
	size_t len = strlen(str);
	if (len == 0) return;
	if (str[len - 1] == '\n') str[len - 1] = '\0';
}
