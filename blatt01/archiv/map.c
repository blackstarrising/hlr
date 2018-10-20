#include <stdio.h>
#include <stdlib.h>
// Definieren Sie ein enum cardd
typedef enum cardd{N = 1, E = 2, S = 4, W = 8} cardd;

// Definieren Sie ein 3x3-Array namens map,das Werte vom Typ cardd enthält
static cardd map[3][3];

// Die Funktion set_dir soll an Position x, y den Wert dir in das Array map eintragen
// Überprüfen Sie x und y um mögliche Arrayüberläufe zu verhindern
// Überprüfen Sie außerdem dir auf Gültigkeit
void set_dir (int x, int y, cardd dir)
{
  //Errorhandling um falsche Eingaben vorzubeugen. Im tatsächlichen Anwendungsfall müsste hier natürlich noch ein Fehler ausgegeben werden etc.
  if(!(x<=2 && x>=0)){return;}
  if(!(y<=2 && y>=0)){return;}
  if(!(dir == N || dir == E ||dir == S ||dir == W || dir == N+W || dir == N+E || dir == S+W || dir == S+E)){return;}

  //An die entsprechende Stelle der Map wird dir geschrieben.
  map[x][y] = dir;
}

// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben
void show_map (void)
{
  //Ein Loop über die gesamte map:
  cardd* ptr = (cardd*)map, * end = ptr + 9;
  int count = 1; //Position auf der Map  mus mitgezählt werden um links- und rechtsbündig zu erzielen (relevant für zweistellige Richtungen)
  static char line[] = "_________"; //Chararray welches eine Zeile der Map darstellt.
  while (ptr != end)
	  {
	    cardd val = *ptr++;
	    //Die FOlgenden Zeilen erreichen das Überprüfen wo in der Karte (links, mitte oder rechts) geschrieben werden soll, ohne eine IF Abfrage zu verwenden.
	    //Je nachdem wird dann links- oder rechtsbündig geschrieben.
	    int pos = ((count-1)%3)*4;
	    int left = pos-(pos/8);
	    int right = pos+1-(pos/8);
	    line[left] = '_';
	    line[right] = '_';

	    //Switch Konstruktion um je nach Eingabe in das Chararray zu schreiben.
	    switch (val)
		{
			case N: line[pos]='N';break;
			case E: line[pos]='E';break;
			case S: line[pos]='S';break;
			case W: line[pos]='W';break;
			case N+W: line[left]='N';line[right]='W';break;
			case N+E: line[left]='N';line[right]='E';break;
			case S+W: line[left]='S';line[right]='W';break;
			case S+E: line[left]='S';line[right]='E';break;
			default: line[pos]='0';break;
		}
	    //Print des Chararrays auf die Konsole:
	    if(count%3 == 0)
	      {
		printf("%s\n",line);
	      }
       	count++;
	}
  printf("%s\n","");
}

int main (void)
{
	// In dieser Funktion darf nichts verändert werden!
	set_dir(0, 1, N);
	set_dir(1, 0, W);
	set_dir(1, 4, W);
	set_dir(1, 2, E);
	set_dir(2, 1, S);

	show_map();

	set_dir(0, 0, N|W);
	set_dir(0, 2, N|E);
	set_dir(0, 2, N|S);
	set_dir(2, 0, S|W);
	set_dir(2, 2, S|E);
	set_dir(2, 2, E|W);
	set_dir(1, 3, N|S|E);
	set_dir(1, 1, N|S|E|W);

	show_map();

	return 0;
}
