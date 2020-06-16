# First run

So, we assume you've got Rivet installed and that the `rivet` command is available in your executable path (i.e. the `$PATH` environment variable). Now it's time to run some events, make some analysis data, and plot the results.


## Generating and analysing

The main interface to Rivet is the `rivet` command. We will demonstrate how to use this to analyse HepMC events from a text file in the "IO_GenEvent" HepMC format.

Firstly, we recommend using a named pipe (or 'FIFO') so that your events
don't create a huge file that takes all your disk space. The idea is that the
generator will push events into what looks like a file, and Rivet will read from it. However, this 'file'
is a direct pipe between the two processes (just like the shell pipe in e.g. `echo foo | less`), so no slow
filesystem access needs to take place, no disk space is occupied, and the system will automatically balance
the data flow between the writing and reading processes. All this is completely transparent to the user:
good old Unix!  Here's how you do it (with a fictional generator command, as an example):
```
> mkfifo fifo.hepmc
> my-generator --num-events=500000 --hepmc-output=fifo.hepmc &
> rivet --analysis=ANALYSIS_NAME fifo.hepmc
```

The backgrounding of the generator process -- i.e. the `&` character -- is important: the generator will wait until the `fifo.hepmc` pipe is being read by Rivet, so unless it is backgrounded you will never get the terminal focus back to run `rivet`!

If you have trouble with this FIFO approach, try writing a few -- 10 or so -- events to a normal file and look at the file (e.g. with emacs, less, gedit or any other text editor). If your generator is writing the events to file properly your HepMC file should consist of a short header, followed by lots of largely numeric lines starting with a letter (E, P, V in particular). If you don't see this, something has gone wrong "upstream" and you should fix it before blaming the FIFO or Rivet ;-)


## Checking the data

By default, Rivet outputs its histograms in the YODA text format. To get a clearer view of the data values in histogram bins, you can use the `yoda2flat` script, e.g. `yoda2flat Rivet.yoda - | less`.

There is also a `yoda2root` command which provides the data as ROOT histogram (or `TGraph`) objects for those of a ROOTy disposition.
