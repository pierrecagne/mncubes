# Makefile for documentation. Requires doxygen.

doc:
	doxygen -u
	doxygen

tex: doc
	make -C doc/latex/

clean:
	rm -rf doc/ Doxyfile.bak *.log

.PHONY : clean
