
# Merging MPI-parallelised Rivet runs on HPC clusters

When analysing events produced on multiple cores there are two common appraches 
one could take to merge the Rivet output into a single output file at the end 
of the run:

## File-based merging

The "brute force" approach is to initialise an `AnalysisHandler` per rank
and have each rank write out a YODA file. At the end of the runs the
various YODA files can be merged using the `rivet-merge` script
(see additional file-based merging documentation [here](merging.md)).

## Merging in memory

Disk space is expensive, however, and it might be more attractive to merge 
the output from each individual rank in memory first, such that only a 
single file needs writing out at the very end.
This can be achieved by using the `AnalysisHandler::writeData(stream, format)` 
method to extract the analysis objects from the `AnalysisHandler` into a 
byte stream for any given rank, which can then be sent across nodes 
and gathered on a single rank. Once collected, it's straightforward to merge 
the different streams using `AnalysisHandlers` together with the 
`AnalysisHandler::readData(stream, format, preload)` method.
Here it is necessary to read in the complete stream by setting the `preload` flag to false
which will load every single analysis object from the stream into the `AnalysisHandler`.
The `AnalysisHandler` ojects can the be merged sequentially, followed by a 
reentrant run on the merged `AnalysisHandler` in order to `finalize()` and 
write out the analysis objects to an output file in the end.

In Python, this could look something like this:

```
from mpi4py import MPI
import rivet, io

def processRank(rank):
  ah = rivet.AnalysisHandler("AH%i" % rank)
  # ... analyse some events ...
  ah.finalize()
  out = io.StringIO()
  ah.writeData(out)
  return out

mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()
mpi_size = mpi_comm.Get_size()

res = processRank(mpi_rank)
res = res.getvalue().encode("utf-8")
res = mpi_comm.gather(res)

if mpi_rank == 0:
  ahmerge = rivet.AnalysisHandler("AHMERGE")
  ahmerge.readData(res[0], preload = False)
  for stream in res[1:]:
    ahtemp = rivet.AnalysisHandler("AHTEMP")
    ahtemp.readData(stream, preload = False)
    ahmerge.merge(ahtemp)
  ahmerge.finalize()
  ahmerge.writeData("mpi_merged_output.yoda.gz")
```
