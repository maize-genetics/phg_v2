!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

# Use Docker to run PHG

# Quick Start

1. Download and install Docker following [instructions](https://docs.docker.com/get-docker/) on the docker website.
2. Pull the PHG docker images from Docker Hub.  It is recommended you pull an image based on a specific tag. Knowing which specific PHG image was used can aid with debugging, and make your pipeline more reproducible.  Change the "0.0.40" tag below to match the PHG tag you wish to pull.
```
docker pull maizegenetics/phg:0.0.40
```

If on a Cornell CBSU machine, use docker1 as below:
```
docker1 pull maizegenetics/phg:0.0.40
```


# Details

[![PHGDocker9pt.png](../images/PHGDocker9pt.png)](../images/PHGDockerDiagramLarge.png)

The PHG is deployed as a Docker image. Docker is an open source platform that wraps applications into an image to be deployed on a host system. Docker containers are instances of a Docker image that may run simultaneously on the same host. See the Docker documentation [here](https://docs.docker.com) for details.

The PHG Docker image is built on the GATK image from [docker hub](https://hub.docker.com). To this image several genomic analysis tools are added, including: VIM, BWA, minimap2, jbwa, htslib, bcftools, SAMTools and TASSEL-5 with the PHG jar. In addition, scripts for loading and processing the PHG are added.

To pull the PHG image from docker hub, run the following commands (replace <version tag> with the tagged version you wish to run):

```
docker pull maizegenetics/phg:<version tag>
```

If on a Cornell CBSU machine, use docker1 as below (Replace <version tag> with the tagged version you wish to run):

```
docker1 pull maizegenetics/phg:<version tag>
```

Access to PHG databases in all cases is through a config file that specifies the host:port, user, password, db name, and the database type (sqlite or postgres). An example of the relevant portions of an SQLite config file is below. NOTE: The user and password parameters are not used for SQLite. The sqlite db file must include a path that is relative to its location within the PHG docker scripts. This will be discussed when each script is discussed.

```
host=localHost
user=sqlite
password=sqlite
DB=/tempFileDir/outputDir/phgTestDB_mapq48.db
DBtype=sqlite
```

An example of a postgres config file is below:
```
host=172.17.0.2:5432\n
user=postgres\n
password=postgres_pwd\n
DB=phgdb\n
DBtype=postgres
```

On CBSU NOTE: If unable to pull a docker image from docker hub on a Cornell CBSU machine, check if you are on an enhanced security machine. cbsumm01, 03, 11 and others are. You can see this if you go to basic biohpc reservations pages, and look at the machines you're able to reserve. They will have a note in the description about enhanced security.


# Kitchen Sink

## *Why PostgreSQL?*
PostgreSQL is an industry standard object-relational database system providing reliability and data integrity.  It provides good performance, supports many concurrent users, and allows for user extensions.

## *Why PostgreSQL in a docker?*
If you do not have access to Postgres on your server, you can download a containerized version of postgres. To do this execute the following comand, replacing <version> with the version of your choice.

```

docker pull postgres:<version>
```
 
 Providing a containerized version of PostgreSQL removes the need for the user to intall and maintain their own postgres image.  The database file created and manipulated in the postgres docker can be persisted to a shared volume for processing by multiple applications.


## *Creating the PHG PostgreSQL docker image*
The postgres docker image may be downloaded from docker hub via the following command.  Replace <version> with the postgres version you wish to run.

```

  docker pull postgres:<version>
```

If running on a Cornell CBSU machine, replace "docker" with "docker1".

If the postgres docker has already been installed and you wish to upgrade, stop and remove any containers associated with the image.  Then remove the existing image using the docker "rmi" command.

## *Connecting to a PostgreSQL docker image*
Connecting to a PHG postgres docker image requires two files:

1. a parameters file
2. a container commands file.

The container commands file reads a user supplied parameters file, then creates a container against the postgres docker image using sourced parameters.  The container is created in detached mode, which allows for "docker exec" commands to be run against it.

Here is an example of a container commands file (container_cmdsPostgres.sh):

```
#!python

#!/usr/bin/env bash
set -u
set -e
set -x

# read in parameters
source $1

# run the docker in detach mode
docker run -i --detach --name $DOCKER_CONTAINER_NAME -e POSTGRES_PASSWORD=$POSTGRES_PASSWORD -v $POSTGRES_DATA_MOUNT:/var/lib/postgresql/data -p $DOCKER_DB_PORT:5432 $DOCKER_IMAGE_NAME

# start the docker
docker start $DOCKER_CONTAINER_NAME

# PHG scripts may now be run using the postgres db with the host ip address and DOCKER_DB_PORT assigned above
```
A sample parameters file, with a description of each parameter, is below.  Among other parameters, the defining of a value for the postgres password helps ensure the security of your data.  Note where each parameter below is used in the container commands file above:

```
#!python

# This file provides parameters for the container_cmdsPostgresDocker.sh file.
# Edit this file to reflect your configuration

# Port that is mapped to the docker postgres 5432 port.  If running on a Cornell CBSU machine,
# the mapped port must be between 8009-8019.  Those are the only exposed ports.
DOCKER_DB_PORT="5433"

# Name for the docker container running against the postgres image
CONTAINER_NAME="phg_postgres_container"

# When running on Cornell CBSU machines, a prefix is appended to any docker container created.
# This prefix consists of the user id followed by "__biohpc_", e.g. "lcj34__biohpc-"
# This prefix can remain blank if not running on Cornell CBSU machine
CONTAINER_PREFIX=""
DOCKER_CONTAINER_NAME=${CONTAINER_PREFIX}${CONTAINER_NAME}

# Name of the docker image against which your container will be run
DOCKER_IMAGE_NAME="postgres:11"

# password for postgres db user.  It is recommended this be changed from the default "postgres"
# The  POSTGRES_PASSWORD must remain consistent when attaching to the same persisted database
# This password should appear in the config file used with other PHG scripts for loading data
# to the phg database.
POSTGRES_PASSWORD="yoursecretpassword"

# Place to mount the data.  Data will be mapped to internal postgres /var/lib/postgresql/data directory
# This data will persist to the volume specified in POSTGRES_DATA_MOUNT.
POSTGRES_DATA_MOUNT="/Users/lcj34/postgres_docker/test_postgres_data/data"
```
To set up the postgres docker, execute the container commands files (container_cmdsPostgres.sh above) from the command line. The db is now ready for connections. Follow the rules below for defining the host/port in the config file if you would like to connect to this postgres instance from the PHG pipeline.

## *Rules of Thumb for Connecting to the PostgreSQL docker image for PHG*
In all cases you must ensure the user name/password you used when creating the database initially remains the same across all containers that connect to the PHG postgres docker image, and in all config files used by scripts that access the database.  

### Connecting from an external machine:

* your config file must have the internet address of the host where the docker container resides
* the port must be the mapped port (value of $DOCKER_DB_PORT) from the "docker run --detach ... command" (see -p = $DOCKER_DB_PORT:5432 in the container and parameters files above)
* example:  If the remote machine's ip address is 128.32.5.269 with a mapped port of 5433, the config file would have a host line of "host=128.32.5.269:5433"

### Connecting from a PHG docker script on the same host:

* your config file should have host ip addr of 0.0.0.0:<db_mapped_port>
* if your machine is using a Docker virtual network, the port would instead be the regular postgres port of 5432
* CAVEAT:  This does not always work.  If you are having problems, try using the external machine address, ie host=128.32.5.269:5433 (replacing the IP address and mapped port with your machine's address and mapped port).

### Connecting from a non-docker program on the same host:

* your config file should have host ip addr of 0.0.0.0:<db_mapped_port>

### Connecting docker-docker on a Cornell CBSU machine (both dockers on the same machine):

* your config file should have host ip addr of 0.0.0.0:5432
* If the above doesn't work, use the full external address with mapped port.  This takes you out and back in the network which shouldn't be required, but it works if nothing else does.
* When a postgres docker run has finished creating and loading a new db, "root" will own all the directories.  In order to run subsequent scripts that mount to this data directory you must "claim" and change permissions on your postgres data directory as below:

```
#!python

> docker1 claim
> chmod -R 777 data  
```


# Troubleshooting

1. Docker can be tricky to get running because it requires root access. If you are having problems, double-check that docker has permissions to access root on your machine.
2. On CBSU NOTE: If unable to pull a docker image from docker hub on a Cornell CBSU machine, check if you are on an enhanced security machine. cbsumm01, 03, 11 and others are. You can see this if you go to basic biohpc reservations pages, and look at the machines you're able to reserve. They will have a note in the description about enhanced security.


[Return to Step 0 pipeline](create_phg_step0_main.md)

[Return to PHG Version 1 instructions](../home_variants_in_gvcf_files.md)

[Return to Wiki Home](../home.md)