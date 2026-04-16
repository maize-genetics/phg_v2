!!! warning "Legacy Documentation - PHG Version 1"

    This section contains documentation for **PHG version 1**, which is
    no longer actively developed. It is preserved here for archival and
    historical reference only. If you are looking to use the Practical
    Haplotype Graph, please refer to the [PHG v2 documentation](../index.md),
    which reflects the current version of the software.

## PHG / BrAPI Setup ##
* git clone https://<userid>@bitbucket.org/bucklerlab/phg_webktor_service.git
* cd phg_webktor_service
* git pull (Do this to get all committed code. This repository should be exactly the master branch. NOT your development environment)
* Edit phg_webktor_service/src/main/resources/application.conf
* Edit phg_webktor_service/docker/config.txt
## Compile PHG / BrAPI Server Code ##
* cd phg_webktor_service
* gradle clean shadowJar
* cp build/libs/phg_webktor_service-1.0-SNAPSHOT-all.jar docker
## Stopping PHG / BrAPI Docker Container ##
* docker1 ps
* docker1 kill <container id>
## Removing PHG / BrAPI Docker Image ##
* docker1 rmi biohpc_tmc46/phg_brapi
## Building PHG / BrAPI Docker Image ##
* cd phg_webktor_service/docker
* docker1 build -t phg_brapi /workdir/tmc46/phg_brapi/phg_webktor_service/docker/
## Starting PHG / BrAPI Docker Container ##
* // Port 80 is external.  Port 8080 is internal
* docker1 run --rm -p 80:8080 biohpc_tmc46/phg_brapi
## Testing PHG / BrAPI Server ##
* curl [http://cbsudc01.biohpc.cornell.edu:80/brapi/v2/serverinfo](http://cbsudc01.biohpc.cornell.edu:80/brapi/v2/serverinfo)
* curl [http://cbsudc01.biohpc.cornell.edu:80/brapi/v2/samples](http://cbsudc01.biohpc.cornell.edu:80/brapi/v2/samples)
* curl [http://cbsudc01.biohpc.cornell.edu:80/brapi/v2/samples/11726](http://cbsudc01.biohpc.cornell.edu:80/brapi/v2/samples/11726)