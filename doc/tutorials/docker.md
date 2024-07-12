# Docker instructions

Docker is a *container* system for providing runnable, pre-built Linux environments -- like virtual machines
but a bit more flexible and composable. We provide various Rivet Docker *images* to allow users to get up
and running with functional Rivet and generator installations, without having to personally confront and
overcome any package-build and system-specific issues.

Obtaining Docker is explained here: [https://docs.docker.com/engine/installation/.](https://docs.docker.com/engine/installation/)

*Note that Docker has historically had a security vulnerability that makes it at least difficult to run
safely on shared computing systems. This is changing via user-space container support, but at the time
of writing (2022) you should probably plan to use the Docker images more as a quick way to try Rivet out,
rather than as a solution for production work.*


## Linux prerequisites
*Note:* on Linux systems, make sure to create a group "docker"
and add yourself to it. Otherwise you will need to run all Docker
commands as sudo.

```
$ sudo groupadd docker
$ sudo usermod -aG docker $USER
```

After that you need to `logout` and log back in for the group change to take
effect. Now you are able to use Docker without sudo.

These steps need only be done once.


## Mac prerequisites
Quote from Docker webpage:

*Docker for Mac requires OS X El Capitan 10.11 or newer macOS release running on
a 2010 or newer Mac, with Intel’s hardware support for MMU virtualization. The
app will run on 10.10.3 Yosemite, but with limited support. Please see What to
know before you install for a full explanation and list of prerequisites.*

Please use the Docker webpage to obtain up-to-date information about possible
Mac install issues.


### Obtaining a Docker image for Rivet

The command is simply:

```
$ docker pull hepstore/rivet:X.Y.Z
```

where X.Y.Z is the latest Rivet version. This will download and store the Docker
container in /var/lib/docker on your Linux or Mac system.

For the rest of this document we will refer to the Rivet/image version number as
"X.Y.Z": replace this with the Rivet version code that you are using, cf. the
`pull` command above.

If you prefer, you can omit the `:X.Y.Z` "tag", in which case the untagged image
retrieved and run will be the latest one available on Docker Hub.


## Running the container interactively

The most basic thing to do is to simply run the container as such:

```
$ docker run -it --rm hepstore/rivet:X.Y.Z
```

(The `-it` flags give you an interactive, terminal container session; the `--rm`
helps the container to clean up after itself. On the latter point, Docker does
like to fill up your drive with saved containers and images: it's a good idea to
periodically call `docker system prune` to clean up and get some disk space
back.)

This gives you full access to anything inside the container.
This is probably a good cross-check to see if your docker permissions
are correct.

To see a list of available Rivet analyses you can now do this:

```
$ docker run -i --rm hepstore/rivet:X.Y.Z rivet --list-analyses
```

You can work with Rivet inside the container exactly as if it were natively
installed on your machine.


## Copying files in and out of the container

You will almost certainly want to copy Rivet output files out of the container
after your run, and maybe to pass in custom analysis files or HepMC event files
for building and processing.

The basic way to copy files in and out of the container, to the host machine, is
to use the `docker copy` command. This will require that you know how to refer
to the open container being used for your interactive session: open another
terminal and run

```
$ docker container ls
CONTAINER ID   IMAGE            COMMAND   CREATED         STATUS         PORTS     NAMES
63c8ae378a0a   hepstore/rivet   "bash"    3 seconds ago   Up 2 seconds             determined_jackson
```

This tells you various things about your container, but we particularly need the
randomly generated container name at the end of the row. This name can be used
to identify the container to be copied into and from:

```
$ docker cp myfile.foo determined_jackson:/work/
$ docker cp determined_jackson:/work/myfile.foo ./
```

This is rather cumbersome, though: wouldn't it be nice to have a better way?
Well, there is one: virtual directories. When you create your container with
`docker run`, add a `-v` flag with a colon-separated pair of directories,
respectively giving paths on the host system and in the container:

```
$ docker run -v $PWD:/host -it --rm hepstore/rivet
```

Now the directory `/host` will exist within your container, and is a view on
the host-system directory from which the run was initiated. You can both read
and write files in this directory, and easily copy files between the container
`/work` area and the `/host`, or use `/host` as your work area if preferred.


## Running MC generators in Docker

Most of the time, we analyse MC event samples not from (large!) pre-generated
files, but on-the-fly either by running a generator and streaming its HepMC output
through a FIFO/interprocess-pipe into Rivet, or by calling Rivet programmatically
from inside the generator. The latter method is more efficient, by avoiding the
conversion to and from the persistent data-format.

This means that the event generator needs to be installed within the container: there
is no obvious way to run an external generator and stream its event output into
Docker for Rivet to analyze. You can generally do this by creating a new "inherited"
Docker image that builds on the `hepstore/rivet` one by installing the generators
and other tools that you want. But fortunately we have already done this for you:
instead of running `hepstore/rivet`, try `hepstore/rivet-pythia` and also `rivet-herwig`,
`rivet-mg5amc` and `rivet-sherpa`!


/*

## DEPRECATED: Running Rivet through docker

In the following, a proposed way of working with rivet through docker is given.

To mount your current directory on the host system into the container and making it the current directory inside the container as well and set the same user and group ids as on the host system, we add this to the command line:

```
 -v $PWD:$PWD -w $PWD -u `id -u $USER`:`id -g`
```


A very efficient way of using rivet through docker is to use aliases as the command line does get quite lengthy.

If you set the following aliases in your shell, you have everything set up to run rivet, compile your own analysis code and make plots:

```
alias rivet='docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:X.Y.Z rivet'
alias rivet-mkanalysis='docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:X.Y.Z rivet-mkanalysis'
alias rivet-buildplugin='docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:X.Y.Z rivet-buildplugin'
alias rivet-mkhtml='docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:X.Y.Z rivet-mkhtml'
alias yodamerge='docker run -i --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:X.Y.Z yodamerge'
```

You might want to put these alias definitions into your `~/.bashrc` for persistence.

*Note:* on SELinux systems, an additional docker run flag `--privileged` is sometimes required for read/write permissions:
```
alias rivet='docker run -i --privileged --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:X.Y.Z rivet'
alias rivet-mkanalysis='docker run -i --privileged --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:X.Y.Z rivet-mkanalysis'
alias rivet-buildplugin='docker run -i --privileged --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:X.Y.Z rivet-buildplugin'
alias rivet-mkhtml='docker run -i --privileged --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:X.Y.Z rivet-mkhtml'
alias yodamerge='docker run -i --privileged --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:X.Y.Z yodamerge'
```

You now can use the following commands from your host system's terminal:

* `rivet`              --- this is the main program
* `rivet-mkanalysis`   --- this creates analysis code templates
* `rivet-build`        --- this compiles your own analysis code
* `rivet-mkhtml`       --- this is used for plotting the output histograms

This should allow you to follow eg. the [firstrun.md](https://gitlab.com/hepcedar/rivet/-/blob/release-3-1-x/doc/tutorials/firstrun.md) tutorial on your host system.

*/
