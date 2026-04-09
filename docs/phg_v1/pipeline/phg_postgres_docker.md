!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../../index.md),
    which reflects the current version of the software.

## *Why PostgreSQL?*
PostgreSQL is an industry standard object-relational database system providing reliability and data integrity.  It provides good performance, supports many concurrent users, and allows for user extensions.

## *Why PostgreSQL in a docker?*
Providing a containerized version of PostgreSQL removes the need for the user to intall and maintain their own postgres image.  The database file created and manipulated in the postgres docker can be persisted to a shared volume for processing by multiple applications.

Having postgres in a docker also ensures the PHG pipeline interacts with a known version of PostgreSQL.

## *Creating the PHG PostgreSQL docker image*
The postgres docker image may be downloaded from docker hub via the following command.  Replace <version> with the postgres version you wish to run.

```
#!python

 > docker pull postgres:<version>
```

If running on a Cornell CBSU machine, replace "docker" with "docker1".

If phg_postgres docker has already been installed and you wish to upgrade, stop and remove any containers associated with the image.  Then remove the existing image using the docker "rmi" command.

## *Connecting to a PHG PostgreSQL docker image*
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
A sample parameters file, with a description of each parameter, is below.  Among other parameters, the defining of   a value for the postgres password helps ensure the security of your data.  Note where each parameter below is used in the container commands file above:

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
To setup the postgres docker execute the container commands files (container_cmdsPostgres.sh above) from the command line.  The db is now ready for connections.  To connect to this postgres instance from the PHG pipeline follow the rules below for defining the host/port in the config file.

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