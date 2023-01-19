This is a project I set for myself to try and write a non-linear time-domain
circuit simulator from scratch in the simplest way possible.

It is very bare bones, only supporting resistors, capacitors, inductors,
diodes, BJTs and sources, but anyone with a bit of linear algebra
experience should be able to extend the modelled components further.

The user-interface is a simple command-line program, that takes a file
describing a circuit, and outputs a .csv file showing the voltages and
currents of selected components and nodes over time. All files ending in
.conf are examples of circuit descriptions.

compile the c files like this:
gcc -lm *.c -o circuitsim

windows executable circuitsim.exe provided, compiled with tcc like this:
tcc *.c -o circuitsim.exe

run like this on windows 
./circuitsim.exe .\astable_multivib.conf
or on linux:
./circuitsim astable_multivib.conf

this produces a file named astable_multivib.conf_results.csv. If you
provided the name of a csv file, like this:
./circuitsim astable_multivib.conf out.csv

then that file name will be used for output.



