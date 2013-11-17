libandor-emu.so: andor-emu.c
	$(CC) -shared andor-emu.c -fPIC -Wall -Wextra -o libandor-emu.so
