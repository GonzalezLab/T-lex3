# T-lex3 Docker instructions (in progress)

By using Docker version of T-lex, you can avoid having to install all required
software and you can have a ready to use version of T-lex3 in less than 30 minutes.

## Requirements

In order to run T-lex3 using Docker you will need:

1. Docker: [windows](https://docs.docker.com/docker-for-windows/install/) / [linux](https://docs.docker.com/install/linux/docker-ce/ubuntu/) / [osx](https://docs.docker.com/docker-for-mac/install/)
2. RepBase for RepeatMasker.

## Build container

We don't provide an already built image because the build process downloads and
uses other software and libraries with copyright and special licenses.

But we are providing the original Dockerfile, so you can build your own copy
using your licenses. But take care that redistributing a container of T-lex3
can be against some software's copyright.

1. First of all, you need to place your RepBase copy in the root of this project,
using the name `RepBase.tar.gz`.

2. To build the image, execute:
    ```
    $ docker build . -t tlex
    ```

    It will take several minutes, as it has to download, compile and install
    several libraries and components.

    During this process, your RepBase library will be included in this image.
    **If you want to update your RepBase library you will need to repeat this
    process.**

## Use the docker container

To use the docker container, you need to know some basics of how Docker works:

- All software and your data will be inside a Docker container, meaning that
only the data you share with the container will be available.

- All T-lex3 process will occur on the container inside the `/data` path, meaning
that you will need to put there a folder which has all the input and data needed
to analyse your data.

To run the example included in this repository **after completing the building phase** you can run:

```
$ docker run -it -v path/to/this/repository/example:/data tlex -O example -A 95 -pairends yes -T /data/TElist_example.txt -M /data/TEannotation_example.txt -G /data/genome_example.fa -R /data/fastq_files
```

The docker execution can be summarized as:

```
$ docker run -it -v /path/to/your/project:/data tlex [ARGUMENTS]
```

## Known issues

- On Windows and OSX, Docker does not run natively and can lead to larger
processing times at some points of the pipeline (however, it's not a big
increase). **On Linux it has native performance.**
