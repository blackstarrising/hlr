#include <stdio.h>

// Definieren Sie ein enum cardd
typedef enum cardd {N = 1, E = 2, S = 4, W = 7} cardd; // enum bekommt Werte deren Werte selbst, als auch die tupel-Summen weiterhin eine bijektive Abbildung auf N sind

// Definieren Sie ein 3x3-Array namens map, das Werte vom Typ cardd enthält
static cardd map[3][3];

// Die Funktion set_dir soll an Position x, y den Wert dir in das Array map eintragen
// Überprüfen Sie x und y um mögliche Arrayüberläufe zu verhindern
// Überprüfen Sie außerdem dir auf Gültigkeit
void set_dir (int x, int y, cardd dir)
{
	// Prüfen ob Indizes korrekt
	if (x > 2 || y > 2 || x < 0 || y < 0)
	{
		//printf('Eingabe ungueltig!');
		return;
	}
	// Prüfen ob cardd korrekt
	if (!(dir == N || dir == E || dir == S || dir == W || dir == N + E || dir == N + W || dir == S + W || dir == S + E))
	{
		//printf("Eingabe ungueltig!");
		return;
	}
	// Setzen von dir in map.
	map[x][y] = dir;
}

// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben
void show_map (void)
{
	for(int y = 0; y <= 2; y++)
	{
		for(int x = 0; x <= 2; x++)
		{
				switch (map[x][y]) {
					case N: printf('%-5c', 'N'); break;
					case E: printf('%-5c', 'E'); break;
					case S: printf('%-5c', 'S'); break;
					case W: printf('%-5c', 'W'); break;
					case N+E: printf('%-1c%-4c', 'N' ,'E'); break;
					case N+W: printf('%-1c%-4c', 'N', 'W'); break;
					case S+E: printf('%-1c%-4c','S', 'E'); break;
					case S+W: printf('%-1c%-4c','S', 'W'); break;
					default: printf('%-5c', '0'); break;
				}
		}
		printf('\n');
	}
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
