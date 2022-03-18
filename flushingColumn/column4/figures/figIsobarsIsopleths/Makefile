 
.SUFFIXES:
.SUFFIXES: .c .cpp .o


UserSolwcad: UserSolwcad.o Solwcad.o NumericalRecipes.o
	g++ -O3 -o UserSolwcad UserSolwcad.o Solwcad.o NumericalRecipes.o 


.cpp.o:
	g++ -O3 -c $<

.c.o:
	gcc -O3 -c $<
