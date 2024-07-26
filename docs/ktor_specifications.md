# Ktor Specifications

In this document, we will discuss the [BrAPI](https://brapicore21.docs.apiary.io/#) 
endpoints provided by the Ktor server implemented in the 
`start-server` command.

## Overview
As stated in the ["Export data"](export_data.md) documentation,
`start-server` is a command that allows users to start a 
"[REST](https://en.wikipedia.org/wiki/REST)ful" web service. This 
service can provide access to a centralized PHG database, allowing 
multiple individuals in a team to simultaneously retrieve 
PHG-relevant information. This service relies on two components:
the Ktor framework and BrAPI-compliant endpoints.

[Ktor](https://ktor.io/) is an asynchronous framework for building 
web applications and microservices in Kotlin. It is designed for 
creating connected systems, focusing on scalability, high 
performance, and ease of use. Ktor supports both server and 
client-side development, allowing developers to create RESTful APIs, 
web sockets, and various other server applications. It provides a 
DSL (Domain-Specific Language) to define routing, request handling, 
and response production, making it highly customizable and flexible 
for a wide range of use cases.

This framework aids in the deployment of BrAPI-compliant data 
endpoints. BrAPI is a standardized RESTful web service API 
specification for plant breeding data. It ensures interoperability 
among breeding databases and tools by providing a common framework. 
Additionally (_and most importantly_), BrAPI is developed by a global 
community of contributors and is intended to be an open and 
accessible standard for anyone involved in plant breeding data 
management.


## Relevant endpoints
While BrAPI comprises several "modules", PHGv2 leverages data
endpoints found within the [core](https://brapicore21.docs.apiary.io/) 
and [genotyping](https://brapigenotyping21.docs.apiary.io/) modules:

* `serverinfo`
    + _Usage_: `<host-url>:<port>/brapi/v2/serverinfo`
    + Find all available BrAPI calls implemented for the server
* `samples`
    + _Usage_: `<host-url>:<port>/brapi/v2/samples`
    + List all samples available in the PHG database
* `variants`
    + _Usage_: `<host-url>:<port>/brapi/v2/variants`
    + List all reference ranges available in the PHG database
* `variantsets`
    + _Usage_: `<host-url>:<port>/brapi/v2/variantsets`
    + Downloadable link for a composite [hVCF](hvcf_specifications.md)
      file. (_Currently for **all** data_)


## Response fields
Each of the prior BrAPI calls will contain "key-value" response
fields which contain the actual data:

* `serverinfo`
    + `calls` - array of available calls to the server
        + `dataTypes` - possible data formats returned by the 
          available call (e.g., `"APPLICATION_JSON"`)
        + `methods` - HTTP methods to be used for each call (e.g., `GET`)
        + `service` - name of the call (e.g., `samples`)
        + `version` - supported versions for a given call (e.g., `"_2"`)
* `samples`
    + `additionalInfo` - "free" space for further information (not
      bound by any data type)
    + `sampleDbId` - internal ID value for given sample
    + `sampleDescription` - description of sample
    + `sampleName` - name of sample in PHG database
* `variants`
    + `additionalInfo` - an object for additional information
    + `alternateBases` - an array of alternate base sequences
    + `ciend` - an array of confidence interval end positions
    + `cipos` - an array of confidence interval start positions
    + `created` - the creation timestamp of the variant
    + `end` - the end position of the variant
    + `filtersApplied` - a boolean indicating if filters were applied
    + `filtersFailed` - an array of failed filters
    + `filtersPassed` - a boolean indicating if the variant passed filters
    + `referenceBases` - the reference base sequence
    + `referenceName` - the reference chromosome or scaffold name
    + `start` - the start position of the variant
    + `svlen` - the structural variant length
    + `updated` - the last update timestamp of the variant
    + `variantDbId` - the unique identifier of the variant
    + `variantNames` - an array of variant names
    + `variantSetDbId` - an array of associated variant set IDs
    + `variantType` - the type of variant (e.g., `"REF_RANGE"`)
* `variantsets`
    + `additionalInfo` - an object for additional information
    + `analysis` - an array of analysis related to the variant set
    + `availableFormats` - a list of available formats with details
        + `dataFormat` - the data format (e.g., `"VCF"`)
        + `fileFormat` - the file format (e.g., `"TEXT_TSV"`)
        + `fileURL` - the URL to access the file
    + `callSetCount` - the count of call sets
    + `referenceSetDbId` - the unique identifier of the reference set
    + `studyDbId` - the unique identifier of the study
    + `variantCount` - the count of variants
    + `variantSetDbId` - the unique identifier of the variant set
    + `variantSetName` - the name of the variant set





