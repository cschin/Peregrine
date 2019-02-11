gcc -g -shared -fPIC sketch.c kalloc.c  -o sketch.so
gcc -g -shared -fPIC mm_select.c  kalloc.c -o mm_select.so
