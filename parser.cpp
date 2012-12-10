/*
 * QDC (quick direct-method controlled) is an optimized exact
 * implementation of the Gillespie's direct-method. It is designed
 * for biochemical simulations when there is the need of dynamic
 * parameters whose values can change during the simulation.
 * version 1.3.4
 *
 * Copyright (C) 2009-2012 Claudio Felicioli
 * mail: c.felicioli@1d20.net - pangon@gmail.com
 *
 * QDC is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

//format note; all the strings needs to be free of characters: ' ' '\t' '\n' ',' '>' '+'

#define MAX_VARIABLES 1000
#define MAX_VARIABLE_NAME_SIZE 100
#define MAX_RULES 1000
#define RULE_TYPE_INITVALUE 0
#define RULE_TYPE_TEMPORALASSIGNMENT 1
#define RULE_NOFPARAMETERS 3
#define MAX_RULES_SIZE 100
#define MAX_REAGENTS 1000
#define MAX_REAGENT_NAME_SIZE 100
#define MAX_REACTIONS 1000
#define MAX_REACTION_FROM 20
#define MAX_REACTION_TO 20
#define MAX_INTRODUCTIONS 1000
#define MAX_TERMNINATION_CONDITIONS 1000

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <limits.h>

typedef long long int int64;
typedef long double double96;

int64 nOfVariables=0;
char variablesNames[MAX_VARIABLES][MAX_VARIABLE_NAME_SIZE];
double96 startingVariablesValues[MAX_VARIABLES];
int64 nOfRules=0;
double96 rules[MAX_RULES][1+RULE_NOFPARAMETERS];
int64 nOfReagents=0;
char reagentsNames[MAX_REAGENTS][MAX_REAGENT_NAME_SIZE];
int64 startingReagents[MAX_REAGENTS];
int64 reagentsIntroduction[MAX_INTRODUCTIONS];
double96 reagentsIntroductionWhen[MAX_INTRODUCTIONS];
int64 reagentsIntroductionHowmany[MAX_INTRODUCTIONS];
int64 nOfIntroductions=0;
int64 nOfReactions=0;
double96 reactK[MAX_REACTIONS];
char reactK_tmp[MAX_REACTIONS][MAX_VARIABLE_NAME_SIZE];
int64 reactFrom[MAX_REACTIONS][MAX_REACTION_FROM];
int64 reactTo[MAX_REACTIONS][MAX_REACTION_TO];
int64 nOfInstReactions=0;
int64 instReactFrom[MAX_REACTIONS][MAX_REACTION_FROM];
int64 instReactFromN[MAX_REACTIONS][MAX_REACTION_FROM];
int64 instReactTo[MAX_REACTIONS][MAX_REACTION_TO];
int64 instReactToN[MAX_REACTIONS][MAX_REACTION_TO];
char* baseline;
char* line;
double96 volume;
double96 time_limit;
bool debug_engine;
int64 nOfMax_values=0;
int64 max_values_reagent[MAX_TERMNINATION_CONDITIONS];
int64 max_values_count[MAX_TERMNINATION_CONDITIONS];
int64 nOfMin_values=0;
int64 min_values_reagent[MAX_TERMNINATION_CONDITIONS];
int64 min_values_count[MAX_TERMNINATION_CONDITIONS];

void trim(char*& s) {
	while(s[0]==' ' || s[0]=='\t' || s[0]=='\n') s++;
	while(strlen(s)>1 && (s[strlen(s)-1]==' ' || s[strlen(s)-1]=='\t' || s[strlen(s)-1]=='\n')) {
		s[strlen(s)-1]='\0';
		}
	}

int64 reagentNumber(char* name) {
	if(strcasecmp("NULL", name)==0) return(-1);
	for(int64 i=0;i<nOfReagents;i++) if(strcmp(reagentsNames[i], name)==0) return(i);
	fprintf(stderr, "format error: \"%s\" undefined.\n", name);
	exit(EXIT_FAILURE);
	}

int64 variableNumber(char* name) {
	for(int64 i=0;i<nOfVariables;i++) if(strcmp(variablesNames[i], name)==0) return(i);
	fprintf(stderr, "format error: \"%s\" undefined.\n", name);
	exit(EXIT_FAILURE);
	}

void readline(FILE* input, const char* errormessage) {
	line=baseline;
	if(fgets(line, 9999, input));
	line[9999]='\0';
	if(feof(input)) {
		if(errormessage!=NULL) fprintf(stderr, "wrong format: %s.\n", errormessage);
		fprintf(stderr, "wrong format\n");
		exit(EXIT_FAILURE);
		}
	trim(line);
	if(line[0]=='#') readline(input, errormessage);
	else {
		char* firstcomm=strchr(line, '#');
		if(firstcomm!=NULL) *firstcomm='\0';
		trim(line);
		}
	}

int main(int nOfArgs, char** args) {
	//type precision checks
	if(sizeof(int64)<8) {
		printf("FATAL ERROR: int64 type bytes < 8 (%zu)\n", sizeof(int64));
		exit(EXIT_FAILURE);
		}
	if(sizeof(double96)<12) {
		printf("FATAL ERROR: double96 type bytes < 12 (%zu)\n", sizeof(double96));
		exit(EXIT_FAILURE);
		}

	//default debug mode
	debug_engine=false;

	//input parameters parsing
	bool recognized_parameter;
	for(int i=1;i<nOfArgs-1;i++) {
		recognized_parameter=false;

		if(strcmp(args[i], "-debug")==0) {
			recognized_parameter=true;
			debug_engine=true;
			}

		if(!recognized_parameter) {
			fprintf(stderr, "unexpected parameter: %s\n", args[i]);
			fprintf(stderr, "usage: %s [-debug] <descriptionfile>\n", args[0]);
			exit(EXIT_FAILURE);
			}
		}

	FILE* input=fopen(args[nOfArgs-1], "r");
	if(input==NULL) {
		fprintf(stderr, "cannot open description file \"%s\".\n", args[nOfArgs-1]);
		exit(EXIT_FAILURE);
		}
	baseline=new char[10000];
	char* token;
	char* a;
	char* token2;
	char* aa;
	char* token3;
	char* aaa;
	char* token4;
	char* aaaa;

	printf("maximum numerical value: %lld\n", LLONG_MAX);

	printf("parsing description file\n");

	printf("parsing reagents declaration\n");
	//reagents
	readline(input, "reagents declaration missing");
	token=strtok_r(line, ",", &a);
	while(token!=NULL) {
		trim(token);
		snprintf(reagentsNames[nOfReagents], MAX_REAGENT_NAME_SIZE-1, "%s", token);
		reagentsNames[nOfReagents][MAX_REAGENT_NAME_SIZE-1]='\0';
		nOfReagents++;
		token=strtok_r(NULL, ",", &a);
		}

	//void line
	readline(input, "missing blank line after reagents declaration");

	printf("parsing volume declaration\n");
	//volume
	readline(input, "volume declaration missing");
	token=strtok_r(line, ",", &a);
	if(token==NULL) {
		fprintf(stderr, "missing or wrong volume parameter.\n");
		exit(EXIT_FAILURE);
		}
	trim(token);
	if(strcmp(token, "volume")!=0) {
		fprintf(stderr, "missing or wrong volume parameter.\n");
		exit(EXIT_FAILURE);
		}
	token=strtok_r(NULL, ",", &a);
	trim(token);
	volume=atof(token);

	//void line
	readline(input, "missing blank line after volume declaration");

	printf("parsing termination conditions declaration\n");
	//time + void line
	readline(input, "time declaration missing");
	token=strtok_r(line, ">", &a);
	if(token==NULL) {
		fprintf(stderr, "missing or wrong time parameter.\n");
		exit(EXIT_FAILURE);
		}
	trim(token);
	if(strcmp(token, "time")!=0) {
		fprintf(stderr, "missing or wrong time parameter.\n");
		exit(EXIT_FAILURE);
		}
	token=strtok_r(NULL, ">", &a);
	trim(token);
	time_limit=atof(token);
	readline(input, "termination conditions declaration must be followed by a blank line");
	bool max_value;
	while(strlen(line)>0) {
		max_value=(strchr(line, '>')!=NULL);

		if(max_value) {
			token=strtok_r(line, ">", &a);
			if(token==NULL) {
				fprintf(stderr, "missing or wrong termination condition declaration (\"%s\").\n", line);
				exit(EXIT_FAILURE);
				}
			trim(token);
			max_values_reagent[nOfMax_values]=reagentNumber(token);
			token=strtok_r(NULL, ">", &a);
			trim(token);
			if(strtoll(token, (char**)NULL, 10)==LLONG_MAX) {
				fprintf(stderr, "WARNING: NUMBERS IN MODEL ARE GREATER THAN MAXIMUM %lld.\n", LLONG_MAX);
				}
			max_values_count[nOfMax_values]=strtoll(token, (char**)NULL, 10);
			nOfMax_values++;
			}
		else {
			token=strtok_r(line, "<", &a);
			if(token==NULL) {
				fprintf(stderr, "missing or wrong termination condition declaration (\"%s\").\n", line);
				exit(EXIT_FAILURE);
				}
			trim(token);
			min_values_reagent[nOfMin_values]=reagentNumber(token);
			token=strtok_r(NULL, "<", &a);
			trim(token);
			if(strtoll(token, (char**)NULL, 10)==LLONG_MAX) {
				fprintf(stderr, "WARNING: NUMBERS IN MODEL ARE GREATER THAN MAXIMUM %lld.\n", LLONG_MAX);
				}
			min_values_count[nOfMin_values]=strtoll(token, (char**)NULL, 10);
			nOfMin_values++;
			}

		readline(input, "termination conditions declaration must be followed by a blank line");
		}

	printf("parsing reactions declaration\n");
	//reactions + void line
	int64 reactFromSize;
	int64 reactToSize;
	readline(input, "reactions declaration missing");
	while(strlen(line)>0) {
		token=strtok_r(line, ",", &a);
		trim(token);
		if(token[0]=='-') {
			token=strtok_r(NULL, ",", &a);
			trim(token);
			token2=strtok_r(token, ">", &aa);
			trim(token2);
			reactFromSize=0;
			token3=strtok_r(token2, "+", &aaa);
			while(token3!=NULL) {
				trim(token3);
//				printf("token \"%s\"\n", token3);
				if(token3[0]>='0' && token3[0]<='9') {
					token4=strtok_r(token3, " ", &aaaa);
//					printf("token4 \"%s\"\n", token4);
					trim(token4);
					if(strtoll(token4, (char**)NULL, 10)==LLONG_MAX) {
						fprintf(stderr, "WARNING: NUMBERS IN MODEL ARE GREATER THAN MAXIMUM %lld.\n", LLONG_MAX);
						}
					instReactFromN[nOfInstReactions][reactFromSize]=strtoll(token4, (char**)NULL, 10);
					token4=strtok_r(NULL, " ", &aaaa);
					if(token4==NULL) {
						token4=token3;
						trim(token4);
						token4+=int(floor(log10(instReactFromN[nOfInstReactions][reactFromSize]))+1);
//						printf("token4a \"%s\"\n", token4);
						}
					else {
						trim(token4);
//						printf("token4b \"%s\"\n", token4);
						}
					}
				else {
					instReactFromN[nOfInstReactions][reactFromSize]=1;
					token4=token3;
					}
				instReactFrom[nOfInstReactions][reactFromSize]=reagentNumber(token4);
				if(++reactFromSize>=MAX_REACTION_FROM) {
					fprintf(stderr, "PARSER ERROR: reagents number exceeds the maximum allowed number MAX_REACTION_FROM=%d\nif really needed you can increase the parameter in the source file parser.cpp\n", MAX_REACTION_FROM);
					exit(EXIT_FAILURE);
					}
				token3=strtok_r(NULL, "+", &aaa);
				}
			instReactFrom[nOfInstReactions][reactFromSize]=-1;
			token2=strtok_r(NULL, ">", &aa);
			trim(token2);
			reactToSize=0;
			token3=strtok_r(token2, "+", &aaa);
			while(token3!=NULL) {
				trim(token3);
				if(token3[0]>='0' && token3[0]<='9') {
					token4=strtok_r(token3, " ", &aaaa);
					trim(token4);
					if(strtoll(token4, (char**)NULL, 10)==LLONG_MAX) {
						fprintf(stderr, "WARNING: NUMBERS IN MODEL ARE GREATER THAN MAXIMUM %lld.\n", LLONG_MAX);
						}
					instReactToN[nOfInstReactions][reactToSize]=strtoll(token4, (char**)NULL, 10);
					token4=strtok_r(NULL, " ", &aaaa);
					if(token4==NULL) {
						token4=token3;
						trim(token4);
						token4+=int(floor(log10(instReactToN[nOfInstReactions][reactToSize]))+1);
//						printf("token4a \"%s\"\n", token4);
						}
					else {
						trim(token4);
//						printf("token4b \"%s\"\n", token4);
						}
					}
				else {
					instReactToN[nOfInstReactions][reactToSize]=1;
					token4=token3;
					}
				instReactTo[nOfInstReactions][reactToSize]=reagentNumber(token4);
				if(++reactToSize>=MAX_REACTION_TO){
					fprintf(stderr, "PARSER ERROR: products number exceeds the maximum allowed number MAX_REACTION_TO=%d\nif really needed you can increase the parameter in the source file parser.cpp\n", MAX_REACTION_TO);
					exit(EXIT_FAILURE);
					}
				token3=strtok_r(NULL, "+", &aaa);
				}
			instReactTo[nOfInstReactions][reactToSize]=-1;
			nOfInstReactions++;
			}
		else {
			if(token[0]=='$') {
//				reactK[nOfReactions]=-1.0-variableNumber(token);
				snprintf(reactK_tmp[nOfReactions], MAX_VARIABLE_NAME_SIZE-1, "%s", token);
				}
			else {
				reactK[nOfReactions]=atof(token);
				reactK_tmp[nOfReactions][0]='#';
				}
			token=strtok_r(NULL, ",", &a);
			trim(token);
			token2=strtok_r(token, ">", &aa);
			trim(token2);
			reactFromSize=0;
			token3=strtok_r(token2, "+", &aaa);
			while(token3!=NULL) {
				trim(token3);
				if(token3[0]=='2') {
					token3++;
					trim(token3);
					reactFrom[nOfReactions][reactFromSize]=reagentNumber(token3);
					if(++reactFromSize>=MAX_REACTION_FROM){
						fprintf(stderr, "PARSER ERROR: reagents number exceeds the maximum allowed number MAX_REACTION_FROM=%d\nif really needed you can increase the parameter in the source file parser.cpp\n", MAX_REACTION_FROM);
						exit(EXIT_FAILURE);
						}
					reactFrom[nOfReactions][reactFromSize]=reagentNumber(token3);
					if(++reactFromSize>=MAX_REACTION_FROM){
						fprintf(stderr, "PARSER ERROR: reagents number exceeds the maximum allowed number MAX_REACTION_FROM=%d\nif really needed you can increase the parameter in the source file parser.cpp\n", MAX_REACTION_FROM);
						exit(EXIT_FAILURE);
						}
					}
				else {
					reactFrom[nOfReactions][reactFromSize]=reagentNumber(token3);
					if(++reactFromSize>=MAX_REACTION_FROM){
						fprintf(stderr, "PARSER ERROR: reagents number exceeds the maximum allowed number MAX_REACTION_FROM=%d\nif really needed you can increase the parameter in the source file parser.cpp\n", MAX_REACTION_FROM);
						exit(EXIT_FAILURE);
						}
					}
				token3=strtok_r(NULL, "+", &aaa);
				}
			reactFrom[nOfReactions][reactFromSize]=-1;
			token2=strtok_r(NULL, ">", &aa);
			trim(token2);
			reactToSize=0;
			token3=strtok_r(token2, "+", &aaa);
			int tmp_number;
			while(token3!=NULL) {
				trim(token3);
				if(token3[0]>='0' && token3[0]<='9') {
					token4=strtok_r(token3, " ", &aaaa);
					trim(token4);
					if(strtoll(token4, (char**)NULL, 10)==LLONG_MAX) {
						fprintf(stderr, "WARNING: NUMBERS IN MODEL ARE GREATER THAN MAXIMUM %lld.\n", LLONG_MAX);
						}
					tmp_number=strtoll(token4, (char**)NULL, 10);
					token4=strtok_r(NULL, " ", &aaaa);
					if(token4==NULL) {
						token4=token3;
						trim(token4);
						token4+=int(floor(log10(tmp_number))+1);
//						printf("token4a \"%s\"\n", token4);
						}
					else {
						trim(token4);
//						printf("token4b \"%s\"\n", token4);
						}
					}
				else {
					tmp_number=1;
					token4=token3;
					}
				for(int tmp_i=0;tmp_i<tmp_number;tmp_i++) {
					reactTo[nOfReactions][reactToSize]=reagentNumber(token4);
					if(++reactToSize>=MAX_REACTION_TO){
						fprintf(stderr, "PARSER ERROR: products number exceeds the maximum allowed number MAX_REACTION_TO=%d\nif really needed you can increase the parameter in the source file parser.cpp\n", MAX_REACTION_TO);
						exit(EXIT_FAILURE);
						}
					}
				token3=strtok_r(NULL, "+", &aaa);
				}

			reactTo[nOfReactions][reactToSize]=-1;
			nOfReactions++;
			}
		readline(input, "reactions declaration must be followed by a blank line");
		}

	//starting reagents + void line
	for(int64 i=0;i<nOfReagents;i++) startingReagents[i]=0;
//	int targetToSet;
	printf("parsing reagent introductions declaration\n");
	readline(input, "missing reagent introductions declaration");
/*	while(strlen(line)>0) {
		token=strtok_r(line, ",", &a);
		trim(token);
		targetToSet=reagentNumber(token);
		token=strtok_r(NULL, ",", &a);
		trim(token);
		startingReagents[targetToSet]=strtoll(token, (char**)NULL, 10);
		readline(input);
		}*/
	while(strlen(line)>0) {
		token=strtok_r(line, ",", &a);
		trim(token);
		reagentsIntroduction[nOfIntroductions]=reagentNumber(token);
		token=strtok_r(NULL, ",", &a);
		trim(token);
		reagentsIntroductionWhen[nOfIntroductions]=atof(token);
		token=strtok_r(NULL, ",", &a);
		trim(token);
		if(strtoll(token, (char**)NULL, 10)==LLONG_MAX) {
			fprintf(stderr, "WARNING: NUMBERS IN MODEL ARE GREATER THAN MAXIMUM %lld.\n", LLONG_MAX);
			}
		reagentsIntroductionHowmany[nOfIntroductions]=strtoll(token, (char**)NULL, 10);
		if(reagentsIntroductionWhen[nOfIntroductions]==0.0) {
			startingReagents[reagentsIntroduction[nOfIntroductions]]=reagentsIntroductionHowmany[nOfIntroductions];
			}
		else nOfIntroductions++;
		readline(input, "missing tailing blank line");
		}

	bool varblock=true;
	REREAD:
	line=baseline;
	if(fgets(line, 9999, input));
	line[9999]='\0';
	if(feof(input)) {
		line[0]='\0';
		varblock=false;
		}
	else {
		trim(line);
		if(line[0]=='#') goto REREAD;
		else {
			char* firstcomm=strchr(line, '#');
			if(firstcomm!=NULL) *firstcomm='\0';
			trim(line);
			}
		}

	printf("parsing variables declaration\n");
	if(varblock) {
		//variables
		//readline(input);
		token=strtok_r(line, ",", &a);
		while(token!=NULL) {
			trim(token);
			if(token[0]!='$') {
				fprintf(stderr, "wrong variable format.\n");
				exit(EXIT_FAILURE);
				}
			snprintf(variablesNames[nOfVariables], MAX_VARIABLE_NAME_SIZE-1, "%s", token);
			variablesNames[nOfVariables][MAX_VARIABLE_NAME_SIZE-1]='\0';
			nOfVariables++;
			token=strtok_r(NULL, ",", &a);
			}
	
		//void line
		readline(input, "missing blank line after variables declaration");
	
		printf("parsing rules declaration\n");
		//rules + void line
		for(int64 i=0;i<nOfVariables;i++) startingVariablesValues[i]=0.0;
		readline(input, "rules declaration missing");
		while(strlen(line)>0) {
			token=strtok_r(line, ",", &a);
			trim(token);
			int64 varid;
			if(strtoll(token, (char**)NULL, 10)==LLONG_MAX) {
				fprintf(stderr, "WARNING: NUMBERS IN MODEL ARE GREATER THAN MAXIMUM %lld.\n", LLONG_MAX);
				}
			switch(strtoll(token, (char**)NULL, 10)) {
				case RULE_TYPE_INITVALUE:
					token=strtok_r(NULL, ",", &a);
					trim(token);
					varid=variableNumber(token);
					token=strtok_r(NULL, ",", &a);
					trim(token);
					startingVariablesValues[varid]=atof(token);
					break;
				case RULE_TYPE_TEMPORALASSIGNMENT:
					rules[nOfRules][0]=RULE_TYPE_TEMPORALASSIGNMENT;
					token=strtok_r(NULL, ",", &a);
					trim(token);
					rules[nOfRules][1]=atof(token);
					token=strtok_r(NULL, ",", &a);
					trim(token);
					rules[nOfRules][2]=variableNumber(token);
					token=strtok_r(NULL, ",", &a);
					trim(token);
					rules[nOfRules][3]=atof(token);
					nOfRules++;
					break;
				default:
					fprintf(stderr, "unknown rule type %lld.\n", strtoll(token, (char**)NULL, 10));
					exit(EXIT_FAILURE);
				}
			readline(input, "missing the blank line after rules declaration");
			}
		}

	for(int64 i=0;i<nOfReactions;i++) {
		if(reactK_tmp[i][0]!='#') {
			reactK[i]=-1.0-variableNumber(reactK_tmp[i]);
			}
		}

	//writing logs to terminal
	for(int64 i=0;i<nOfReagents;i++) printf("reagent %lld=%s : %lld\n", i, reagentsNames[i], startingReagents[i]);
	printf("volume=%.300Lf\n", volume);
	for(int64 i=0;i<nOfVariables;i++) printf("variable %lld=%s : %Lf\n", i, variablesNames[i], startingVariablesValues[i]);
	for(int64 i=0;i<nOfRules;i++) {
		switch(int64(rules[i][0])) {
			case RULE_TYPE_INITVALUE:
				fprintf(stderr, "parse error.\n");
				exit(EXIT_FAILURE);
				break;
			case RULE_TYPE_TEMPORALASSIGNMENT:
				printf("at time %Lf variable %s=%Lf\n", rules[i][1], variablesNames[int64(rules[i][2])], rules[i][3]);
				break;
			default:
				fprintf(stderr, "unknown rule type %lld.\n", int64(rules[i][0]));
				exit(EXIT_FAILURE);
			}
		}
	for(int64 i=0;i<nOfReactions;i++) {
		if(reactK[i]<0) {
			printf("reaction %lld (k=%s) : ", i, variablesNames[int64(-reactK[i])-1]);
			}
		else printf("reaction %lld (k=%Lf) : ", i, reactK[i]);
		int64 j=0;
		if(reactFrom[i][j]==-1) {
			printf("NULL");
			}
		else {
			while(reactFrom[i][j]!=-1) {
				if(j>0) printf(" + ");
				printf("%s", reagentsNames[reactFrom[i][j]]);
				j++;
				}
			}
		printf(" > ");
		j=0;
		if(reactTo[i][j]==-1) {
			printf("NULL");
			}
		else {
			while(reactTo[i][j]!=-1) {
				if(j>0) printf(" + ");
				printf("%s", reagentsNames[reactTo[i][j]]);
				j++;
				}
			}
		printf("\n");
		}

	//generating src/reactions.cpp
	printf("generating \"src/reactions.cpp\".\n");

	FILE* output=fopen("src/reactions.cpp", "w");
	if(output==NULL) {
		fprintf(stderr, "error: fopen the output file fails\n");
		exit(EXIT_FAILURE);
		}

	fprintf(output, "/*\n * QDC (quick direct-method controlled) is an optimized exact\n * implementation of the \
Gillespie's direct-method. It is designed\n * for biochemical simulations when there is the need of dynamic\n * \
parameters whose values can change during the simulation.\n * version 1.3.4\n *\n * Copyright (C) 2009-2012 Claudio \
Felicioli\n * mail: c.felicioli@1d20.net - pangon@gmail.com\n *\n * QDC is free software; you can redistribute it \
and/or modify\n * it under the terms of the GNU General Public License as published by\n * the Free Software \
Foundation; either version 3 of the License, or\n * (at your option) any later version.\n *\n * This program is \
distributed in the hope that it will be useful,\n * but WITHOUT ANY WARRANTY; without even the implied warranty \
of\n * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n * GNU General Public License for more details.\
\n *\n * You should have received a copy of the GNU General Public License\n * along with this program. If not, \
see <http://www.gnu.org/licenses/>.\n */\n\n");

	for(int64 i=0;i<nOfReactions;i++) {
		//detecting reaction order
		int64 k=0;
		while(reactFrom[i][k]!=-1) k++;
		if(k==0) {
			if(reactK[i]<0) {
				fprintf(output, "#define basalRate%lld Avo*globalvar_%s\n", i, variablesNames[int64(-reactK[i])-1]+1);
				}
			else fprintf(output, "#define basalRate%lld Avo*%.16Lf\n", i, reactK[i]);
			}
		else if(k==1) {
			if(reactK[i]<0) {
				fprintf(output, "#define basalRate%lld globalvar_%s\n", i, variablesNames[int64(-reactK[i])-1]+1);
				}
			else fprintf(output, "#define basalRate%lld %.16Lf\n", i, reactK[i]);
			}
		else if(k==2) {
			if(reactK[i]<0) {
				fprintf(output, "#define basalRate%lld (globalvar_%s/(Vol*Avo))\n", i, variablesNames[int64(-reactK[i])-1]+1);
				}
			else fprintf(output, "#define basalRate%lld (%.16Lf/(Vol*Avo))\n", i, reactK[i]);
			}
		else {
			fprintf(stderr, "error: unimplemented %lld-order\n", k);
			exit(EXIT_FAILURE);
			}
		//generating reaction execution function
		fprintf(output, "void reaction%lld(void) {\n", i);
		int64 j=0;
		while(reactFrom[i][j]!=-1) {
			fprintf(output, "\treagent%lld(-1);\n", reactFrom[i][j]);
			j++;
			}
		j=0;
		while(reactTo[i][j]!=-1) {
			fprintf(output, "\treagent%lld(1);\n", reactTo[i][j]);
			j++;
			}
		fprintf(output, "\t}\n\n");
		}

	for(int64 i=0;i<nOfInstReactions;i++) {
		fprintf(output, "void instReaction%lld(void) {\n", i);
		int64 j=0;
		while(instReactFrom[i][j]!=-1) {
			fprintf(output, "\treagent%lld(-%lldLL);\n", instReactFrom[i][j], instReactFromN[i][j]);
			j++;
			}
		j=0;
		while(instReactTo[i][j]!=-1) {
			fprintf(output, "\treagent%lld(%lldLL);\n", instReactTo[i][j], instReactToN[i][j]);
			j++;
			}
		fprintf(output, "\t}\n\n");
		}

	//closing file
	int output_fd;
	if((output_fd=fileno(output))==-1) {
		fprintf(stderr, "failed getting output file descriptor");
		exit(EXIT_FAILURE);
		}
	if(fsync(output_fd)!=0) {
		fprintf(stderr, "failed fsync of output");
		exit(EXIT_FAILURE);
		}
	if(fclose(output)!=0) {
		fprintf(stderr, "failed fclose of output");
		exit(EXIT_FAILURE);
		}

	//generating src/reagents.cpp
	printf("generating \"src/reagents.cpp\".\n");

	output=fopen("src/reagents.cpp", "w");
	if(output==NULL) {
		fprintf(stderr, "error: fopen the output file fails\n");
		exit(EXIT_FAILURE);
		}

	fprintf(output, "/*\n * QDC (quick direct-method controlled) is an optimized exact\n * implementation of the \
Gillespie's direct-method. It is designed\n * for biochemical simulations when there is the need of dynamic\n * \
parameters whose values can change during the simulation.\n * version 1.3.4\n *\n * Copyright (C) 2009-2012 Claudio \
Felicioli\n * mail: c.felicioli@1d20.net - pangon@gmail.com\n *\n * QDC is free software; you can redistribute it \
and/or modify\n * it under the terms of the GNU General Public License as published by\n * the Free Software \
Foundation; either version 3 of the License, or\n * (at your option) any later version.\n *\n * This program is \
distributed in the hope that it will be useful,\n * but WITHOUT ANY WARRANTY; without even the implied warranty \
of\n * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n * GNU General Public License for more details.\
\n *\n * You should have received a copy of the GNU General Public License\n * along with this program. If not, \
see <http://www.gnu.org/licenses/>.\n */\n\n");

	bool inTheLeft;
	bool twoInTheLeft;
	for(int64 i=0;i<nOfReagents;i++) {
		//generating reagent modification function
		fprintf(output, "void reagent%lld(int64 mod) {\n", i);
		fprintf(output, "\treagent[%lld]+=mod;\n", i);
		//calculating propensity modifications
		for(int64 j=0;j<nOfReactions;j++) {
			inTheLeft=false;
			twoInTheLeft=false;
			int64 k=0;
			while(reactFrom[j][k]!=-1) {
				if(reactFrom[j][k]==i) {
					if(inTheLeft) twoInTheLeft=true;
					inTheLeft=true;
					}
				k++;
				}
			if(inTheLeft) {
				if(k==1) fprintf(output, "\tau[%lld]=reagent[%lld]*basalRate%lld;\n", j, i, j);
				else if(k==2) {
					if(twoInTheLeft)fprintf(output, "\tau[%lld]=(basalRate%lld*reagent[%lld]*(reagent[%lld]-1))/2;\n", j, j, reactFrom[j][0], reactFrom[j][1]);
					else fprintf(output, "\tau[%lld]=basalRate%lld*reagent[%lld]*reagent[%lld];\n", j, j, reactFrom[j][0], reactFrom[j][1]);
					}
				else {
					fprintf(stderr, "error: unimplemented %lld-order\n", k);
					exit(EXIT_FAILURE);
					}
				if(debug_engine) {
					fprintf(output, "\t//DEBUG\n");
					fprintf(output, "\tif(au[%lld]<0.0) {\n", j);
					fprintf(output, "\t\tprintf(\"error: a0=%%Lf\\n\", a0);\n");
					fprintf(output, "\t\tfor(int64 j=0;j<nOfReagents;j++) printf(\"%%lld,\", reagent[j]);\n");
					fprintf(output, "\t\tprintf(\"\\n\");\n");
					fprintf(output, "\t\tfor(int64 j=0;j<nOfReactions;j++) printf(\"%%Lf,\", au[j]);\n");
					fprintf(output, "\t\tprintf(\"\\n\");\n");
					fprintf(output, "\t\texit(EXIT_FAILURE);\n");
					fprintf(output, "\t\t}\n");
					}
				}
			}
		for(int64 j=0;j<nOfInstReactions;j++) {
			inTheLeft=false;
			int64 k=0;
			while(instReactFrom[j][k]!=-1) {
				if(instReactFrom[j][k]==i) inTheLeft=true;
				k++;
				}
			if(inTheLeft) {
				fprintf(output, "\tinst[%lld]=", j);
				k=0;
				while(instReactFrom[j][k]!=-1) {
					if(k>0) fprintf(output, " && ");
					fprintf(output, "(reagent[%lld] >= %lldLL)", instReactFrom[j][k], instReactFromN[j][k]);
					k++;
					}
				fprintf(output, ";\n");
				}
			}
		fprintf(output, "\t}\n\n");
		}

	//closing file
	if((output_fd=fileno(output))==-1) {
		fprintf(stderr, "failed getting output file descriptor");
		exit(EXIT_FAILURE);
		}
	if(fsync(output_fd)!=0) {
		fprintf(stderr, "failed fsync of output");
		exit(EXIT_FAILURE);
		}
	if(fclose(output)!=0) {
		fprintf(stderr, "failed fclose of output");
		exit(EXIT_FAILURE);
		}

	//generating src/init.cpp
	printf("generating \"src/init.cpp\".\n");

	output=fopen("src/init.cpp", "w");
	if(output==NULL) {
		fprintf(stderr, "error: fopen the output file fails\n");
		exit(EXIT_FAILURE);
		}

	fprintf(output, "/*\n * QDC (quick direct-method controlled) is an optimized exact\n * implementation of the \
Gillespie's direct-method. It is designed\n * for biochemical simulations when there is the need of dynamic\n * \
parameters whose values can change during the simulation.\n * version 1.3.4\n *\n * Copyright (C) 2009-2012 Claudio \
Felicioli\n * mail: c.felicioli@1d20.net - pangon@gmail.com\n *\n * QDC is free software; you can redistribute it \
and/or modify\n * it under the terms of the GNU General Public License as published by\n * the Free Software \
Foundation; either version 3 of the License, or\n * (at your option) any later version.\n *\n * This program is \
distributed in the hope that it will be useful,\n * but WITHOUT ANY WARRANTY; without even the implied warranty \
of\n * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n * GNU General Public License for more details.\
\n *\n * You should have received a copy of the GNU General Public License\n * along with this program. If not, \
see <http://www.gnu.org/licenses/>.\n */\n\n");

	fprintf(output, "void init(void) {\n");
	fprintf(output, "\tprintf(\"control int64 size   : %%d bytes\\n\", sizeof(int64));\n");
	fprintf(output, "\tprintf(\"control double96 size: %%d bytes\\n\", sizeof(double96));\n");
	fprintf(output, "\tif(sizeof(int64)<8 || sizeof(double96)<12) {\n");
	fprintf(output, "\t\tfprintf(stderr, \"FATAL ERROR: wrong type sizes\\n\");\n");
	fprintf(output, "\t\texit(EXIT_FAILURE);\n");
	fprintf(output, "\t\t}\n");
	fprintf(output, "\tdouble v=Vol;\n");
	fprintf(output, "\tprintf(\"verify of volume [double]  : %%.300f\\n\", v);\n");
	fprintf(output, "\tdouble96 v96=Vol;\n");
	fprintf(output, "\tprintf(\"verify of volume [double96]: %%.300Lf\\n\", v96);\n");
	for(int64 i=0;i<nOfReactions;i++) {
		fprintf(output, "\treactionDo[%lld]=(void*)(&reaction%lld);\n", i, i);
		}
	for(int64 i=0;i<nOfInstReactions;i++) {
		fprintf(output, "\tinstReactionDo[%lld]=(void*)(&instReaction%lld);\n", i, i);
		}
	fprintf(output, "\tfor(int64 i=0;i<nOfReactions;i++) au[i]=0.0;\n");
	fprintf(output, "\tfor(int64 i=0;i<nOfInstReactions;i++) inst[i]=false;\n");
	fprintf(output, "\tfor(int64 i=0;i<nOfReagents;i++) reagent[i]=0;\n");
	for(int64 i=0;i<nOfReagents;i++) fprintf(output, "\treagent%lld(%lldLL);\n", i, startingReagents[i]);
	for(int64 i=0;i<nOfReactions;i++) {
		if(reactFrom[i][0]==-1) fprintf(output, "\tau[%lld]=basalRate%lld;\n", i, i);
		}
	fprintf(output, "\t}\n\n");

	//closing file
	if((output_fd=fileno(output))==-1) {
		fprintf(stderr, "failed getting output file descriptor");
		exit(EXIT_FAILURE);
		}
	if(fsync(output_fd)!=0) {
		fprintf(stderr, "failed fsync of output");
		exit(EXIT_FAILURE);
		}
	if(fclose(output)!=0) {
		fprintf(stderr, "failed fclose of output");
		exit(EXIT_FAILURE);
		}

	//generating src/engine.cpp
	printf("generating \"src/engine.cpp\".\n");

	output=fopen("src/engine.cpp", "w");
	if(output==NULL) {
		fprintf(stderr, "error: fopen the output file fails\n");
		exit(EXIT_FAILURE);
		}
	fprintf(output, "/*\n * QDC (quick direct-method controlled) is an optimized exact\n * implementation of the \
Gillespie's direct-method. It is designed\n * for biochemical simulations when there is the need of dynamic\n * \
parameters whose values can change during the simulation.\n * version 1.3.4\n *\n * Copyright (C) 2009-2012 Claudio \
Felicioli\n * mail: c.felicioli@1d20.net - pangon@gmail.com\n *\n * QDC is free software; you can redistribute it \
and/or modify\n * it under the terms of the GNU General Public License as published by\n * the Free Software \
Foundation; either version 3 of the License, or\n * (at your option) any later version.\n *\n * This program is \
distributed in the hope that it will be useful,\n * but WITHOUT ANY WARRANTY; without even the implied warranty \
of\n * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n * GNU General Public License for more details.\
\n *\n * You should have received a copy of the GNU General Public License\n * along with this program. If not, \
see <http://www.gnu.org/licenses/>.\n */\n\n");

	fprintf(output, "#define nOfReagents %lld\n", nOfReagents);
	fprintf(output, "#define nOfReactions %lld\n", nOfReactions);
	fprintf(output, "#define nOfInstReactions %lld\n", nOfInstReactions);
//	fprintf(output, "#define Vol 0.0000000000000000003\n\n");
	fprintf(output, "#define Vol %.300LfL\n", volume);
	fprintf(output, "#define Avo 602200000000000000000000.0L\n\n");

	fprintf(output, "typedef long long int int64;\n");
	fprintf(output, "typedef long double double96;\n\n");

	fprintf(output, "#include <stdlib.h>\n#include <stdio.h>\n#include <time.h>\n#include <math.h>\n#include <unistd.h>\n#include <string.h>\n\n");
	for(int64 i=0;i<nOfVariables;i++) fprintf(output, "double96 globalvar_%s=%.16Lf;\n", variablesNames[i]+1, startingVariablesValues[i]);
	for(int64 i=0;i<nOfReagents;i++) fprintf(output, "void reagent%lld(int64 mod);\n", i);
	fprintf(output, "void* reactionDo[nOfReactions];\n");
	fprintf(output, "void* instReactionDo[nOfInstReactions];\n");
	fprintf(output, "int64 reagent[nOfReagents];\n");
	fprintf(output, "double96 au[nOfReactions];\n");
	fprintf(output, "bool inst[nOfInstReactions];\n");
	fprintf(output, "int64 instCount;\n");
	fprintf(output, "int64 instSel;\n");
	fprintf(output, "FILE* output;\n");
	fprintf(output, "FILE* output2;\n");
	fprintf(output, "FILE* output3;\n");
	fprintf(output, "FILE* output4;\n");
	fprintf(output, "double96 t;\n");
	fprintf(output, "double96 a;\n");
	fprintf(output, "double96 a0;\n");
	fprintf(output, "double96 at;\n");
	fprintf(output, "double96 tao_denom=0.0;\n");
	fprintf(output, "double96 tao=0.0;\n");
	fprintf(output, "int64 selectedReaction;\n");
	fprintf(output, "int64 reactCounts[nOfReactions+nOfInstReactions];\n");
	fprintf(output, "#include \"reactions.cpp\"\n#include \"reagents.cpp\"\n#include \"init.cpp\"\n\n");
	fprintf(output, "\n");
	fprintf(output, "void GLIteration(void);\n\n");
	fprintf(output, "int main(int nOfArgs, char** args) {\n");
	fprintf(output, "\tif(nOfArgs!=3 && nOfArgs!=2 && nOfArgs!=1) {\n");
	fprintf(output, "\t\tfprintf(stderr, \"usage: %%s [samplingFrequency [randomSeed]]\\n\", args[0]);\n");
	fprintf(output, "\t\texit(EXIT_FAILURE);\n");
	fprintf(output, "\t\t}\n\n");
	fprintf(output, "\tdouble96 tLimit=%.20LfL;\n", time_limit);
	fprintf(output, "\tdouble96 tSampling;\n\n");
	fprintf(output, "\tif(nOfArgs>=2) tSampling=atof(args[1]);\n");
	fprintf(output, "\telse tSampling=0.1;\n\n");
	fprintf(output, "\tfor(int64 i=0;i<(nOfReactions+nOfInstReactions);i++) reactCounts[i]=0;\n");
	fprintf(output, "\tinit();\n\n");
	fprintf(output, "\tif((output=fopen(\"%s_reagents.csv\", \"w\"))==NULL) {\n", args[nOfArgs-1]);
	fprintf(output, "\t\tfprintf(stderr, \"error: fopen the reagents output file fails\\n\");\n");
	fprintf(output, "\t\texit(EXIT_FAILURE);\n");
	fprintf(output, "\t\t}\n");
	fprintf(output, "\tif((output2=fopen(\"%s_reactions.csv\", \"w\"))==NULL) {\n", args[nOfArgs-1]);
	fprintf(output, "\t\tfprintf(stderr, \"error: fopen the reactions output file fails\\n\");\n");
	fprintf(output, "\t\texit(EXIT_FAILURE);\n");
	fprintf(output, "\t\t}\n");
	fprintf(output, "\tif((output3=fopen(\"%s_reactioncounts.csv\", \"w\"))==NULL) {\n", args[nOfArgs-1]);
	fprintf(output, "\t\tfprintf(stderr, \"error: fopen the reactioncounts output file fails\\n\");\n");
	fprintf(output, "\t\texit(EXIT_FAILURE);\n");
	fprintf(output, "\t\t}\n");
	fprintf(output, "\tif((output4=fopen(\"%s_log.txt\", \"w\"))==NULL) {\n", args[nOfArgs-1]);
	fprintf(output, "\t\tfprintf(stderr, \"error: fopen the log output file fails\\n\");\n");
	fprintf(output, "\t\texit(EXIT_FAILURE);\n");
	fprintf(output, "\t\t}\n\n");
	fprintf(output, "\tfprintf(output,\"");
	for(int64 i=0;i<nOfReagents;i++) {
		if(i>0) fprintf(output, ",");
		fprintf(output, "%s", reagentsNames[i]);
		}
	fprintf(output, "\\n\");\n");
	fprintf(output, "\tfprintf(output2,\"");
	for(int64 i=0;i<nOfReactions;i++) {
		if(i>0) fprintf(output, ",");
		int64 j=0;
		if(reactFrom[i][j]==-1) {
			fprintf(output, "NULL");
			}
		else {
			while(reactFrom[i][j]!=-1) {
				if(j>0) fprintf(output, "+");
				fprintf(output, "%s", reagentsNames[reactFrom[i][j]]);
				j++;
				}
			}
		fprintf(output, " > ");
		j=0;
		if(reactTo[i][j]==-1) {
			fprintf(output, "NULL");
			}
		else {
			while(reactTo[i][j]!=-1) {
				if(j>0) fprintf(output, "+");
				fprintf(output, "%s", reagentsNames[reactTo[i][j]]);
				j++;
				}
			}
		}
	fprintf(output, ",a0\\n\");\n");
	fprintf(output, "\tfprintf(output3,\"");
	for(int64 i=0;i<nOfReactions;i++) {
		if(i>0) fprintf(output, ",");
		int64 j=0;
		if(reactFrom[i][j]==-1) {
			fprintf(output, "NULL");
			}
		else {
			while(reactFrom[i][j]!=-1) {
				if(j>0) fprintf(output, "+");
				fprintf(output, "%s", reagentsNames[reactFrom[i][j]]);
				j++;
				}
			}
		fprintf(output, " > ");
		j=0;
		if(reactTo[i][j]==-1) {
			fprintf(output, "NULL");
			}
		else {
			while(reactTo[i][j]!=-1) {
				if(j>0) fprintf(output, "+");
				fprintf(output, "%s", reagentsNames[reactTo[i][j]]);
				j++;
				}
			}
		}
	for(int64 i=0;i<nOfInstReactions;i++) {
		if(i>0 || nOfReactions>0) fprintf(output, ",");
		int64 j=0;
		while(instReactFrom[i][j]!=-1) {
			if(j>0) fprintf(output, "+");
			fprintf(output, "%lld %s", instReactFromN[i][j], reagentsNames[instReactFrom[i][j]]);
			j++;
			}
		fprintf(output, " > ");
		j=0;
		while(instReactTo[i][j]!=-1) {
			if(j>0) fprintf(output, "+");
			fprintf(output, "%lld %s", instReactToN[i][j], reagentsNames[instReactTo[i][j]]);
			j++;
			}
		}
	fprintf(output, "\\n\");\n");

	fprintf(output, "\tfprintf(output4,\"model parsing command parameters:\\n\");\n");
	fprintf(output, "\tfprintf(output4,\"");
	for(int i=0;i<nOfArgs;i++) fprintf(output, "%s ", args[i]);
	fprintf(output, "\\n\\n\");\n");

	fprintf(output, "\tfprintf(output4,\"simulation execution parameters:\\n\");\n");
	fprintf(output, "\tfor(int i=0;i<nOfArgs;i++) fprintf(output4, \"%%s \", args[i]);\n");
	fprintf(output, "\tfprintf(output4,\"\\n\\n\");\n");

	fprintf(output, "\tfprintf(output4,\"model name: %s\\n\\n\");\n", args[nOfArgs-1]);

	fprintf(output, "\ttime_t timer_start;\n");
	fprintf(output, "\tstruct tm* tm_info;\n");
	fprintf(output, "\tchar tm_string[50];\n");
	fprintf(output, "\ttime(&timer_start);\n");
	fprintf(output, "\ttm_info=localtime(&timer_start);\n");
	fprintf(output, "\tstrftime(tm_string, 50, \"%%H:%%M:%%S %%Y/%%m/%%d\", tm_info);\n");
	fprintf(output, "\tfprintf(output4,\"simulation start time: %%s\\n\\n\", tm_string);\n");

	fprintf(output, "\tint randomseed;\n");
	fprintf(output, "\tif(nOfArgs==3) randomseed=atoi(args[2]);\n");
	fprintf(output, "\telse randomseed=time(NULL);\n");
	fprintf(output, "\tsrandom(randomseed);\n");
//	fprintf(output, "\tdouble96 tLimit=40.0;\n");
	//fprintf(output, "\tdouble96 tSampling=0.01;\n");
	fprintf(output, "\tdouble96 tSamplingTemp=0.0;\n");
	fprintf(output, "\tint64 ite=0;\n");
	fprintf(output, "\t//int64 iteLimit=10000000;\n");
	fprintf(output, "\tt=0.0;\n");
	fprintf(output, "\tchar lastLine[10000];\n");
	fprintf(output, "\tchar lastLine2[10000];\n");
	fprintf(output, "\tchar lastLine3[100000];\n\n");
	fprintf(output, "\tlastLine[0]='\\0';\n");
	fprintf(output, "\tfor(int64 i=0;i<nOfReagents;i++) snprintf(lastLine+strlen(lastLine), 9999-strlen(lastLine), \"%%lld,\", reagent[i]);\n");
	fprintf(output, "\tsnprintf(lastLine+strlen(lastLine), 9999-strlen(lastLine), \"\\n\");\n");
	fprintf(output, "\tfprintf(output, \"%%s\", lastLine);\n");
	fprintf(output, "\tlastLine2[0]='\\0';\n");
	fprintf(output, "\tfor(int64 i=0;i<nOfReactions;i++) snprintf(lastLine2+strlen(lastLine2), 9999-strlen(lastLine2), \"%%Lf,\", au[i]);\n");
	fprintf(output, "\ta0=0.0;\n");
	fprintf(output, "\tfor(int64 i=0;i<nOfReactions;i++) a0+=au[i];\n");
	fprintf(output, "\tsnprintf(lastLine2+strlen(lastLine2), 9999-strlen(lastLine2), \"%%Lf,\", a0);\n");
	fprintf(output, "\tsnprintf(lastLine2+strlen(lastLine2), 9999-strlen(lastLine2), \"\\n\");\n");
	fprintf(output, "\tfprintf(output2, \"%%s\", lastLine2);\n");
	fprintf(output, "\tlastLine3[0]='\\0';\n");
	fprintf(output, "\tfor(int64 i=0;i<(nOfReactions+nOfInstReactions);i++) snprintf(lastLine3+strlen(lastLine3), 99999-strlen(lastLine3), \"%%lld,\", reactCounts[i]);\n");
	fprintf(output, "\tsnprintf(lastLine3+strlen(lastLine3), 99999-strlen(lastLine3), \"\\n\");\n");
	fprintf(output, "\tfprintf(output3, \"%%s\", lastLine3);\n\n");
	//fprintf(output, "\twhile(t<tLimit && ite<iteLimit) {\n");
	fprintf(output, "\twhile(t<tLimit) {\n");	
	for(int64 i=0;i<nOfRules;i++) {
		switch(int64(rules[i][0])) {
			case RULE_TYPE_INITVALUE:
				fprintf(stderr, "parse error.\n");
				exit(EXIT_FAILURE);
				break;
			case RULE_TYPE_TEMPORALASSIGNMENT:
				fprintf(output, "\t\tif((t-tao)<%.16Lf && t>=%.16Lf) {\n\t\t\tglobalvar_%s=%.16Lf;\n", rules[i][1], rules[i][1], variablesNames[int64(rules[i][2])]+1, rules[i][3]);
//				fprintf(output, "\t\tif(t<%.16Lf && t+tao>=%.16Lf) {\n\t\t\tglobalvar_%s=%.16Lf;\n", rules[i][1], rules[i][1], variablesNames[int64(rules[i][2])]+1, rules[i][3]);
				for(int64 j=0;j<nOfReagents;j++) fprintf(output, "\t\t\treagent%lld(1);reagent%lld(-1);\n", j, j);
				fprintf(output, "\t\t\t}\n");
				break;
			default:
				fprintf(stderr, "unknown rule type %lld.\n", int64(rules[i][0]));
				exit(EXIT_FAILURE);
			}
		}

	for(int64 i=0;i<nOfIntroductions;i++) {
		//no perche' potrei volerlo azzerare ad un tempo n, piuttosto sarebbe meglio non specificare nel file descrizione le inizializzazioni a 0
//		if(reagentsIntroductionHowmany[i]!=0) {
			fprintf(output, "\t\tif((t-tao)<%.16Lf && t>=%.16Lf) {\n", reagentsIntroductionWhen[i], reagentsIntroductionWhen[i]);
//			fprintf(output, "\t\tif(t<%.16Lf && t+tao>=%.16Lf) {\n", reagentsIntroductionWhen[i], reagentsIntroductionWhen[i]);
			fprintf(output, "\t\t\treagent%lld(%lldLL);\n", reagentsIntroduction[i], reagentsIntroductionHowmany[i]);
			fprintf(output, "\t\t\t}\n");
//			}
		}

	fprintf(output, "\t\tGLIteration();\n");
	fprintf(output, "\t\tt+=tao;\n");
	fprintf(output, "\t\ttSamplingTemp+=tao;\n");
	fprintf(output, "\t\tite++;\n");
	fprintf(output, "\t\twhile(tSamplingTemp>(tSampling*2.0) && t<tLimit) {\n");
	fprintf(output, "\t\t\tfprintf(output, \"%%s\", lastLine);\n");
	fprintf(output, "\t\t\tfprintf(output2, \"%%s\", lastLine2);\n");
	fprintf(output, "\t\t\tfprintf(output3, \"%%s\", lastLine3);\n");
	fprintf(output, "\t\t\ttSamplingTemp-=tSampling;\n");
	fprintf(output, "\t\t\t}\n");
	fprintf(output, "\t\tif(tSamplingTemp>tSampling) {\n");
	fprintf(output, "\t\t\ttSamplingTemp-=tSampling;\n");
	fprintf(output, "\t\t\tlastLine[0]='\\0';\n");
	fprintf(output, "\t\t\tfor(int64 i=0;i<nOfReagents;i++) snprintf(lastLine+strlen(lastLine), 9999-strlen(lastLine), \"%%lld,\", reagent[i]);\n");
	fprintf(output, "\t\t\tsnprintf(lastLine+strlen(lastLine), 9999-strlen(lastLine), \"\\n\");\n");
	fprintf(output, "\t\t\tfprintf(output, \"%%s\", lastLine);\n");
	fprintf(output, "\t\t\tlastLine2[0]='\\0';\n");
	fprintf(output, "\t\t\tfor(int64 i=0;i<nOfReactions;i++) snprintf(lastLine2+strlen(lastLine2), 9999-strlen(lastLine2), \"%%Lf,\", au[i]);\n");
	fprintf(output, "\t\t\ta0=0.0;\n");
	fprintf(output, "\t\t\tfor(int64 i=0;i<nOfReactions;i++) a0+=au[i];\n");
	fprintf(output, "\t\t\tsnprintf(lastLine2+strlen(lastLine2), 9999-strlen(lastLine2), \"%%Lf,\", a0);\n");
	fprintf(output, "\t\t\tsnprintf(lastLine2+strlen(lastLine2), 9999-strlen(lastLine2), \"\\n\");\n");
	fprintf(output, "\t\t\tfprintf(output2, \"%%s\", lastLine2);\n");
	fprintf(output, "\t\t\tlastLine3[0]='\\0';\n");
	fprintf(output, "\t\t\tfor(int64 i=0;i<(nOfReactions+nOfInstReactions);i++) snprintf(lastLine3+strlen(lastLine3), 99999-strlen(lastLine3), \"%%lld,\", reactCounts[i]);\n");
	fprintf(output, "\t\t\tsnprintf(lastLine3+strlen(lastLine3), 99999-strlen(lastLine3), \"\\n\");\n");
	fprintf(output, "\t\t\tfprintf(output3, \"%%s\", lastLine3);\n");
	fprintf(output, "\t\t\t}\n");
	fprintf(output, "\t\tif(ite%%100000==0) printf(\"iteration %%lld   t=%%.20Lf\\n\", ite, t);\n");

	for(int64 i=0;i<nOfMax_values;i++) fprintf(output, "\t\tif(reagent[%lld]>%lldLL) break;\n", max_values_reagent[i], max_values_count[i]);
	for(int64 i=0;i<nOfMin_values;i++) fprintf(output, "\t\tif(reagent[%lld]<%lldLL) break;\n", min_values_reagent[i], min_values_count[i]);

	fprintf(output, "\t\t}\n\n");

	fprintf(output, "\ttime_t timer_end;\n");
	fprintf(output, "\ttime(&timer_end);\n");
	fprintf(output, "\tfprintf(output4,\"run time: %%ld seconds\\n\\n\", timer_end-timer_start);\n");

	fprintf(output, "\tfprintf(output4,\"simulated time: %%Lf seconds\\n\\n\", t);\n");

	fprintf(output, "\tint64 total_reactCounts=0;\n");
	fprintf(output, "\tfor(int i=0;i<nOfReactions+nOfInstReactions;i++) total_reactCounts+=reactCounts[i];\n");
	fprintf(output, "\tfprintf(output4,\"total reactions count: %%lld reactions\\n\\n\", total_reactCounts);\n");

	fprintf(output, "\tfprintf(output4,\"random seed: %%d\\n\\n\", randomseed);\n\n");

	fprintf(output, "\tint output_fd;\n");
	fprintf(output, "\tif((output_fd=fileno(output))==-1) {\n");
	fprintf(output, "\t\tfprintf(stderr, \"failed getting output file descriptor\");\n");
	fprintf(output, "\t\texit(EXIT_FAILURE);\n");
	fprintf(output, "\t\t}\n");
	fprintf(output, "\tif(fsync(output_fd)!=0) {\n");
	fprintf(output, "\t\tfprintf(stderr, \"failed fsync of output\");\n");
	fprintf(output, "\t\texit(EXIT_FAILURE);\n");
	fprintf(output, "\t\t}\n");
	fprintf(output, "\tif(fclose(output)!=0) {\n");
	fprintf(output, "\t\tfprintf(stderr, \"failed fclose of output\");\n");
	fprintf(output, "\t\texit(EXIT_FAILURE);\n");
	fprintf(output, "\t\t}\n");
	fprintf(output, "\tif((output_fd=fileno(output2))==-1) {\n");
	fprintf(output, "\t\tfprintf(stderr, \"failed getting output file descriptor\");\n");
	fprintf(output, "\t\texit(EXIT_FAILURE);\n");
	fprintf(output, "\t\t}\n");
	fprintf(output, "\tif(fsync(output_fd)!=0) {\n");
	fprintf(output, "\t\tfprintf(stderr, \"failed fsync of output\");\n");
	fprintf(output, "\t\texit(EXIT_FAILURE);\n");
	fprintf(output, "\t\t}\n");
	fprintf(output, "\tif(fclose(output2)!=0) {\n");
	fprintf(output, "\t\tfprintf(stderr, \"failed fclose of output\");\n");
	fprintf(output, "\t\texit(EXIT_FAILURE);\n");
	fprintf(output, "\t\t}\n");
	fprintf(output, "\tif((output_fd=fileno(output3))==-1) {\n");
	fprintf(output, "\t\tfprintf(stderr, \"failed getting output file descriptor\");\n");
	fprintf(output, "\t\texit(EXIT_FAILURE);\n");
	fprintf(output, "\t\t}\n");
	fprintf(output, "\tif(fsync(output_fd)!=0) {\n");
	fprintf(output, "\t\tfprintf(stderr, \"failed fsync of output\");\n");
	fprintf(output, "\t\texit(EXIT_FAILURE);\n");
	fprintf(output, "\t\t}\n");
	fprintf(output, "\tif(fclose(output3)!=0) {\n");
	fprintf(output, "\t\tfprintf(stderr, \"failed fclose of output\");\n");
	fprintf(output, "\t\texit(EXIT_FAILURE);\n");
	fprintf(output, "\t\t}\n");
	fprintf(output, "\tif((output_fd=fileno(output4))==-1) {\n");
	fprintf(output, "\t\tfprintf(stderr, \"failed getting output file descriptor\");\n");
	fprintf(output, "\t\texit(EXIT_FAILURE);\n");
	fprintf(output, "\t\t}\n");
	fprintf(output, "\tif(fsync(output_fd)!=0) {\n");
	fprintf(output, "\t\tfprintf(stderr, \"failed fsync of output\");\n");
	fprintf(output, "\t\texit(EXIT_FAILURE);\n");
	fprintf(output, "\t\t}\n");
	fprintf(output, "\tif(fclose(output4)!=0) {\n");
	fprintf(output, "\t\tfprintf(stderr, \"failed fclose of output\");\n");
	fprintf(output, "\t\texit(EXIT_FAILURE);\n");
	fprintf(output, "\t\t}\n\n");
	fprintf(output, "\texit(EXIT_SUCCESS);\n");
	fprintf(output, "\t}\n\n");
	fprintf(output, "void GLIteration() {\n");
	fprintf(output, "\ta0=0.0;\n");
	fprintf(output, "\tfor(int64 i=0;i<nOfReactions;i++) {\n");
	fprintf(output, "\t\ta0+=au[i];\n");

	if(debug_engine) {
		fprintf(output, "\t\t//DEBUG\n");
		fprintf(output, "\t\tif(a0<0.0) {\n");
		fprintf(output, "\t\t\tprintf(\"error: a0=%%Lf\\n\", a0);\n");
		fprintf(output, "\t\t\tfor(int64 j=0;j<nOfReagents;j++) printf(\"%%lld,\", reagent[j]);\n");
		fprintf(output, "\t\t\tprintf(\"\\n\");\n");
		fprintf(output, "\t\t\tfor(int64 j=0;j<nOfReactions;j++) printf(\"%%Lf,\", au[j]);\n");
		fprintf(output, "\t\t\tprintf(\"\\n\");\n");
		fprintf(output, "\t\t\texit(EXIT_FAILURE);\n");
		fprintf(output, "\t\t\t}\n");
		}
	fprintf(output, "\t\t}\n");
	fprintf(output, "\tif(a0==0.0) {\n");
	fprintf(output, "\t\tprintf(\"error: a0=0.0\\n\");\n");
	fprintf(output, "\t\texit(EXIT_FAILURE);\n");
	fprintf(output, "\t\t}\n");
	fprintf(output, "\ta=a0*(((double)(random()%%RAND_MAX))/((double)RAND_MAX));\n");
	fprintf(output, "\tat=0.0;\n");
	fprintf(output, "\tselectedReaction=0;\n");
	fprintf(output, "\twhile(a>at+au[selectedReaction]) {\n");
	fprintf(output, "\t\tat+=au[selectedReaction];\n");
	fprintf(output, "\t\tselectedReaction++;\n");
	if(debug_engine) {
		fprintf(output, "\t\t//DEBUG\n");
		fprintf(output, "\t\tif(selectedReaction>=nOfReactions) {\n");
		fprintf(output, "\t\t\tprintf(\"error: selectedReaction==nOfReactions (a=%%Lf, a0=%%Lf)\\n\", a, a0);\n");
		fprintf(output, "\t\t\tfor(int64 j=0;j<nOfReagents;j++) printf(\"%%lld,\", reagent[j]);\n");
		fprintf(output, "\t\t\tprintf(\"\\n\");\n");
		fprintf(output, "\t\t\tfor(int64 j=0;j<nOfReactions;j++) printf(\"%%Lf,\", au[j]);\n");
		fprintf(output, "\t\t\tprintf(\"\\n\");\n");
		fprintf(output, "\t\t\texit(EXIT_FAILURE);\n");
		fprintf(output, "\t\t\t}\n");
		}
	fprintf(output, "\t\t}\n");
	fprintf(output, "\tif(au[selectedReaction]==0.0) {\n");
	fprintf(output, "\t\tprintf(\"WARNING: au[%%lld]==0.0 (a=%%Lf, a0=%%Lf)\\n\", selectedReaction, a, a0);\n");
	fprintf(output, "\t\tfor(int64 j=0;j<nOfReagents;j++) printf(\"%%lld,\", reagent[j]);\n");
	fprintf(output, "\t\tprintf(\"\\n\");\n");
	fprintf(output, "\t\tfor(int64 j=0;j<nOfReactions;j++) printf(\"%%Lf,\", au[j]);\n");
	fprintf(output, "\t\tprintf(\"\\n\");\n");
	fprintf(output, "\t\twhile(au[selectedReaction]==0.0) {\n");
	fprintf(output, "\t\t\tselectedReaction++;\n");
	fprintf(output, "\t\t\tif(selectedReaction>=nOfReactions) {\n");
	fprintf(output, "\t\t\t\tselectedReaction--;\n");
	fprintf(output, "\t\t\t\twhile(au[selectedReaction]==0.0) {\n");
	fprintf(output, "\t\t\t\t\tselectedReaction--;\n");
	fprintf(output, "\t\t\t\t\t}\n");
	fprintf(output, "\t\t\t\t}\n");
	fprintf(output, "\t\t\t}\n");
	fprintf(output, "\t\t}\n");
	fprintf(output, "\ttao_denom=((double)(random()%%RAND_MAX))/((double)RAND_MAX);\n");
	fprintf(output, "\twhile(tao_denom==0.0) {\n");
	fprintf(output, "\t\tprintf(\"WARNING: tao_denom==0.0\\n\");\n");
	fprintf(output, "\t\ttao_denom=((double)(random()%%RAND_MAX))/((double)RAND_MAX);\n");
	fprintf(output, "\t\t}\n");
	fprintf(output, "\ttao=(1.0/a0)*logl(1.0/tao_denom);\n");
	if(debug_engine) {
		fprintf(output, "\t//DEBUG\n");
		fprintf(output, "\tif(tao<0.0) {\n");
		fprintf(output, "\t\tprintf(\"tao<0 :%%Lf\\n\", tao);\n");
		fprintf(output, "\t\tfor(int64 j=0;j<nOfReagents;j++) printf(\"%%lld,\", reagent[j]);\n");
		fprintf(output, "\t\tprintf(\"\\n\");\n");
		fprintf(output, "\t\tfor(int64 j=0;j<nOfReactions;j++) printf(\"%%Lf,\", au[j]);\n");
		fprintf(output, "\t\tprintf(\"\\n\");\n");
		fprintf(output, "\t\texit(EXIT_FAILURE);\n");
		fprintf(output, "\t\t}\n");
		}
	fprintf(output, "\t((void (*)(void))(reactionDo[selectedReaction]))();\n");
	fprintf(output, "\treactCounts[selectedReaction]++;\n");
	fprintf(output, "\twhile(true) {\n");
	fprintf(output, "\t\tinstCount=0;\n");
	fprintf(output, "\t\tfor(int64 i=0;i<nOfInstReactions;i++) if(inst[i]) instCount++;\n");
	fprintf(output, "\t\tif(instCount==0) break;\n");
	fprintf(output, "\t\tif(instCount>1) instSel=random()%%instCount;\n");
	fprintf(output, "\t\telse instSel=0;\n");
	fprintf(output, "\t\tfor(int64 i=0;i<nOfInstReactions;i++) {\n");
	fprintf(output, "\t\t\tif(inst[i]) {\n");
	fprintf(output, "\t\t\t\tif(instSel==0) {\n");
	fprintf(output, "\t\t\t\t\t((void (*)(void))(instReactionDo[i]))();\n");
	fprintf(output, "\t\t\t\t\treactCounts[nOfReactions+i]++;\n");
	fprintf(output, "\t\t\t\t\tbreak;\n");
	fprintf(output, "\t\t\t\t\t}\n");
	fprintf(output, "\t\t\t\telse instSel--;\n");
	fprintf(output, "\t\t\t\t}\n");
	fprintf(output, "\t\t\t}\n");
	fprintf(output, "\t\t}\n");
	fprintf(output, "\t}\n\n");

	//closing file
	if((output_fd=fileno(output))==-1) {
		fprintf(stderr, "failed getting output file descriptor");
		exit(EXIT_FAILURE);
		}
	if(fsync(output_fd)!=0) {
		fprintf(stderr, "failed fsync of output");
		exit(EXIT_FAILURE);
		}
	if(fclose(output)!=0) {
		fprintf(stderr, "failed fclose of output");
		exit(EXIT_FAILURE);
		}

	//closing description file
	int input_fd;
	if((input_fd=fileno(input))==-1) {
		fprintf(stderr, "failed getting input file descriptor");
		exit(EXIT_FAILURE);
		}
	if(fsync(input_fd)!=0) {
		fprintf(stderr, "failed fsync of input");
		exit(EXIT_FAILURE);
		}
	if(fclose(input)!=0) {
		fprintf(stderr, "failed fclose of input");
		exit(EXIT_FAILURE);
		}

	//freeing memory
	delete[] baseline;

	exit(EXIT_SUCCESS);
	}
