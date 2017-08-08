:- module(dtw,[
	       dtw/4
	      ]).

% author: Vito Cuccovillo

:- ensure_loaded(library(lists)).

%---[DTW]--- INIZIO complessita' O(n^2)
% formati in input accettati:
% dtw([1,1,3,1,3,1,1],[1,3,1,2,1,1,1],S,M).
% dtw([a,1,b,3,c,2,d,3,e],[c,2,d,2,e,4,f],S,M).
% questo è il formalismo senza lo zero iniziale! dunque la sequenza finisce con un termine
% dietro le quinte inserisco comunque lo zero iniziale per avere tutte le corrispondenze dtw
% formalismo dei mappings in output: [[a-N,b-M]...] necessario per capire univocamente le corrispondenze
% a è l'elemento della prima seq, N è l'indice di tale elemento della seq
% b è l'elemento nella seconda seq, M è l'indice dell'elemento nella seconda seq
%print('WARPING MATRIX'),nl,print_matrix(M5),
dtw([],[],0,[]) :- !. % caso di path vuoti, restituisce vuoto
dtw(A,[],0,[]) :- !. %idem
dtw([],B,0,[]) :- !.  %idem
dtw(Lx,Ly,S,M) :-
	extract_weights(Lx,L11),
	extract_weights(Ly,L22),
	(isSequenceWithLabels(Lx) -> % se ho nomi di eventi aggiungo lo 0 iniziale
		appendZero(L11,L1), % necessario per il dtw
		appendZero(L22,L2) % necessario per il dtw
		;
		L1 = L11,
		L2 = L22
	),
	length(L1,A), 
	length(L2,B), 
	L1L is A+1,
	L2L is B+1, 
	init_dtw(L1L,L2L,M1),
	firstColInf(M1,L1L,M2), 
	firstRowInf(M2,L2L,M3), 
	replace(M3,0,0,0,M4), 
	process(M4,1,1,L1,L2,M5),	
	getPath(M5,A,B,L1,L2,M6), 
	at(M5,A,B,L), 
	length(M6,LP), 
	S is L/LP,
	reverse(M6,M7), % riordina gli accoppiamenti
	extract_labels(Lx,Labelx),
	extract_labels(Ly,Labely),
	get_mapping_labels(M7,Labelx,Labely,M).

isSequenceWithLabels([X,L|T]) :-
	not(number(X)),
	number(L).

% da una lista di eventi con peso [a,1,b,3,c,1] 
% inserisco come peso del primo termine lo zero, per segnare l'inizio (utile per dtw)
% restituisce una lista di soli pesi in base alle posizioni
extract_weights([],[]).
extract_weights([X],[]) :- % caso termine-peso, con termine finale
	not(number(X)).
extract_weights([X,N|T],[N|T1]) :- % caso termine-peso
	number(N),
	not(number(X)),
	!,
	extract_weights(T,T1).
extract_weights([N|T],[N|T1]) :- % caso solo peso
	number(N),
	!,
	extract_weights(T,T1).


% da una lista di eventi con peso
% restituisce una lista di sole etichette rispettando le posizioni
extract_labels([],[]).
extract_labels([X],[X]).
extract_labels([X,N|T],[X|T1]) :- % caso termine-peso
	number(N),
	not(number(X)),
	!,
	extract_labels(T,T1).
extract_labels([N|T],[N|T1]) :- % caso solo peso
	!,
	extract_labels(T,T1).

appendZero(L,[0|L]).

% restituisce gli accoppiamenti del DTW indicando le etichette
% invece che le posizioni degli elementi nella lista
% output: [etichetta1-indice1, etichetta2-indice2]
get_mapping_labels([],Lx,Ly,[]).
get_mapping_labels([[X,Y]|T],Lx,Ly,[[Ex-X,Ey-Y]|T1]) :-
	atV(X,Lx,Ex), 
	atV(Y,Ly,Ey),
	!,
	get_mapping_labels(T,Lx,Ly,T1).

process(M,I,J,L1,L2,N) :- 
	length(L2,A), 
	B is A+1, 
	J < B, 
	II is I -1, JJ is J -1, 
	atV(II,L1,X),
	atV(JJ,L2,Y),
	sequence_distance(X,Y,Cost),
	offset(M,I,J,O), 
	NewCost is Cost + O,
	replace(M,I,J,NewCost,N1),
	length(L2,L2L), 
	NJ is J + 1, 
	!,
	process(N1,I,NJ,L1,L2,N).

process(M,I,J,L1,L2,N) :- 
	length(L1,A), 
	B is A+1, 
	I < B,
	II is I+1, 
	!,
	process(M,II,1,L1,L2,N).
process(M,I,J,L1,L2,M).

offset(M,I,J,O) :-
	IM is I-1, 
	JM is J-1,
	at(M,IM,J,M1),
	at(M,I,JM,M2),
	at(M,IM,JM,M3),
	minimum(M1,M2,M3,O).

%restituisce gli abbinamenti degli oggetti più simili percorrendo a ritroso il miglior percorso nella matrice
getPath(M,A,B,L1,L2,[]) :- 
	A = 0, B = 0.
getPath(M,X,Y,L1,L2,[[II,JJ]|T]) :- 
	II is X-1, JJ is Y-1,
	at(M,X,Y,O), 
	nextStep(M,X,Y,NX,NY),
	!,
	getPath(M,NX,NY,L1,L2,T).	

nextStep(M,X,Y,NX,NY) :- 
	MX is X-1, 
	MY is Y-1,
	at(M,MX,MY,C1), 
	at(M,X,MY,C2), 
	at(M,MX,Y,C3),
	minPath([[MX,MY,C1],[X,MY,C2],[MX,Y,C3]],NX,NY).

% minimo fra 3 valori per DTW
minimum(A,B,C,M) :- 
	P is min(A,B), 
	M is min(P,C).

% trova la casella migliore andando a ritroso
minPath([[X1,Y1,C1],[X2,Y2,C2],[X3,Y3,C3]],NX,NY) :-
	C1 < C2, 
	!,
	minPath([[X1,Y1,C1],[X3,Y3,C3]],NX,NY).
minPath([[X1,Y1,C1],[X2,Y2,C2],[X3,Y3,C3]],NX,NY) :-
	!,
	minPath([[X2,Y2,C2],[X3,Y3,C3]],NX,NY).
minPath([[X2,Y2,C2],[X3,Y3,C3]],X2,Y2) :-
	C2 < C3.
minPath([[X2,Y2,C2],[X3,Y3,C3]],X3,Y3).
minPath([[X1,Y1,C1]],X1,Y1).

% creare L1 righe con L2 valori (inizializzati a sup (+inf)).
init_dtw(0,L2,[]).
init_dtw(L1,L2,[L|T]) :- 
	N is L1 - 1, 
	createRow(L2,L),
	!,
	init_dtw(N,L2,T).

%metti ad infinito la prima colonna '+Inf'
firstColInf(L,1,N) :- 
	replace(L,0,0,inf,N).
firstColInf(M,L,N) :-
	I is L -1, 
	replace(M,I,0,inf,N1), 
	!,
	firstColInf(N1,I,N).

%metti ad infinito la prima colonna
firstRowInf(L,1,N) :- 
	replace(L,0,0,inf,N).
firstRowInf(M,L,N) :-
	I is L -1, 
	replace(M,0,I,inf,N1), 
	!,
	firstRowInf(N1,I,N).

% crea una riga con N volte 0
createRow(0,[]).
createRow(N,[0|T]) :- 
	M is N - 1, 
	!,
	createRow(M,T).

sequence_distance(A,B,D) :-
	N is A-B, 
	abs(N,M), 
	D is M*M.

abs(X,X) :- 
	X >= 0.
abs(X,Y) :- 
	X < 0, 
	Y is -X.
	
printList([]).
printList([A|B]) :-
  format('~w\t',A),
  printList(B).

print_matrix([]) :- 
	nl.
print_matrix([H|T]) :- 
	write(H), 
	nl, 
	print_matrix(T).

%---[DTW]--- FINE

%PREDICATI GESTIONE MATRICI

row(M, N, Row) :-
    nth1(N, M, Row).

column(M, N, Col) :-
    transpose(M, MT),
    row(MT, N, Col).

symmetrical(M) :-
    transpose(M, M).

transpose([[]|_], []) :- !.
transpose([[I|Is]|Rs], [Col|MT]) :-
    first_column([[I|Is]|Rs], Col, [Is|NRs]),
	!,
    transpose([Is|NRs], MT).

first_column([], [], []).
first_column([[]|_], [], []).
first_column([[I|Is]|Rs], [I|Col], [Is|Rest]) :-
	!,
    first_column(Rs, Col, Rest).
	
% valore specificando le coordinate	
% at lavora 1-based, il sistema invece 0-based, dunque aggiungo 1 alle coordinate
at(Mat, Row, Col, Val) :- 
	R is Row+1, 
	C is Col+1,
	nth1(R, Mat, ARow), 
	nth1(C, ARow, Val).
	
atV(N,L,X) :- 
	nth0(N,L,X).

%Replace: sostituisce il valore alla posizione (X,Y) con il nuovo valore
% e restituisce la nuova matrice con la sostituzione
replace([L|Ls],0,Y,Z,[R|Ls]) :- % della riga indicata,
  replace_column(L,Y,Z,R),
  !.                % si sostituisce la colonna
  
replace([L|Ls],X,Y,Z,[L|Rs]) :- % se non trovi la riga
	X > 0 ,                                 % e l'indice è positivo
	X1 is X-1 ,								% decrementa l'indice
	!, 
	replace(Ls,X1,Y,Z,Rs). 
                           
replace_column([_|Cs],0,Z,[Z|Cs]) :- !.  %trovato l'offset fai la sostituzione
replace_column([C|Cs],Y,Z,[C|Rs]) :- % altrimenti
	Y > 0 ,                                    			% se l'indice è positivo
	Y1 is Y-1 ,  										% decrementalo
	!,  
	replace_column(Cs,Y1,Z,Rs).         

	