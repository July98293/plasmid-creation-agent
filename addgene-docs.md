openapi: 3.0.3
info:
title: Addgene Developers API
version: e88b4172e8
description: |2

    ### Welcome to the Addgene Developers API

    This documentation describes the Addgene Developers API and how to interact with it programmatically.
    The examples throughout use [curl](https://curl.se/download.html), and we recommend [Bruno](https://docs.usebruno.com/bruno-basics/download) as a graphical client for exploring and testing endpoints.

    You can download or inspect the raw OpenAPI schema using the **Download** button at the top of this page, or access it directly at: https://api.developers.addgene.org/docs/schema/.

    ### Tools

    - **curl** — Command-line examples are written using curl. You can download it from [curl.se](https://curl.se/download.html); it is often preinstalled on macOS and Linux.
    - **Bruno** — A lightweight graphical API client. After installing [Bruno](https://docs.usebruno.com/bruno-basics/download), create a new collection and import the OpenAPI schema from https://api.developers.addgene.org/docs/schema/.
    Bruno makes it easy to send requests, inspect responses, and manage authentication tokens in one place.

paths:
/auth/test/:
get:
operationId: auth_test_retrieve
description: |2

        ### Summary

        Endpoint to test authentication.

        Requires scope `auth:test`.

        ### Example Usage

        ```bash
        curl --request GET \
          --header "Authorization: Token ${KEY}" \
          --url https://api.developers.addgene.org/auth/test/
        ```
      tags:
      - auth
      security:
      - tokenAuth: []
      responses:
        '200':
          content:
            application/json:
              schema:
                $ref: '#/components/schemas/TestAuthResponse'
              examples:
                SuccessfulResponse:
                  value:
                    username: testuser
                    email: testuser@example.com
                    auth_method: Token
                  summary: Example of a successful response
                  description: This example shows the structure of a successful response.
          description: ''

/download/plasmids/:
get:
operationId: download_plasmids_retrieve
description: |2+

        ### Summary

        Endpoint to download plasmid data.

        Requires scope `bulk-download:plasmids`.

        ### Example Usage

        ```bash
        curl --request GET \
          --location \
          --header "Authorization: Token ${KEY}" \
          --url https://api.developers.addgene.org/download/plasmids/
        ```


        ### File Structure

        The data is provided as JSON. JSON has the following properties:

        - It is encoded with UTF-8.
        - It isn't ordered, so the keys can appear in any order.
        - Whitespace/line breaks may vary.
        - Any array can contain zero or more items.
        - Fields inapplicable to a given item may be null.
        - All keys are included, even if they hold empty string or null values.

        JSON parsers and libraries exist for almost every programming language. We recommend you use one.

        *\* Principal Investigator/Head of Laboratory names are not included in the data.*


        ### File Versioning

        Our JSON file follows a version format: `majorVersion.minorVersion`.

        Changes in the `majorVersion` may require you to update your code, while `minorVersion` updates are generally backward-compatible.
        Specialized versions may have a suffix format (e.g., `1.0-XXX`), while the main feed will have no suffix.

      tags:
      - download
      security:
      - tokenAuth: []
      responses:
        '200':
          content:
            application/json:
              schema:
                $ref: '#/components/schemas/PlasmidsList'
          description: ''

/download/plasmids_with_sequences/:
get:
operationId: download_plasmids_with_sequences_retrieve
description: |2+

        ### Summary

        Endpoint to download plasmids with sequences data.

        Requires scope `bulk-download:plasmids-with-sequences`.

        ### Example Usage

        ```bash
        curl --request GET \
          --location \
          --header "Authorization: Token ${KEY}" \
          --url https://api.developers.addgene.org/download/plasmids_with_sequences/
        ```


        ### File Structure

        The data is provided as JSON. JSON has the following properties:

        - It is encoded with UTF-8.
        - It isn't ordered, so the keys can appear in any order.
        - Whitespace/line breaks may vary.
        - Any array can contain zero or more items.
        - Fields inapplicable to a given item may be null.
        - All keys are included, even if they hold empty string or null values.

        JSON parsers and libraries exist for almost every programming language. We recommend you use one.

        *\* Principal Investigator/Head of Laboratory names are not included in the data.*


        ### File Versioning

        Our JSON file follows a version format: `majorVersion.minorVersion`.

        Changes in the `majorVersion` may require you to update your code, while `minorVersion` updates are generally backward-compatible.
        Specialized versions may have a suffix format (e.g., `1.0-XXX`), while the main feed will have no suffix.

      tags:
      - download
      security:
      - tokenAuth: []
      responses:
        '200':
          content:
            application/json:
              schema:
                $ref: '#/components/schemas/PlasmidsWithSequencesList'
          description: ''

/download/genbank/{id}/:
get:
operationId: download_genbank_retrieve
description: |2

        ### Summary

        Endpoint to download a GenBank file for a given sequence id.

        Requires scope `bulk-download:plasmids-with-sequences` or `catalog:retrieve-with-sequences`.

        ### Example Usage

        ```bash
        curl --request GET \
          --location \
          --header "Authorization: Token ${KEY}" \
          --url https://api.developers.addgene.org/download/genbank/438456/
        ```
      parameters:
      - in: path
        name: id
        schema:
          type: integer
        description: Sequence ID
        required: true
      tags:
      - download
      security:
      - tokenAuth: []
      responses:
        '302':
          description: No response body

/catalog/plasmid/:
get:
operationId: catalog_plasmid_list
description: |2

        ### Summary

        Endpoint to retrieve a list of plasmids.

        Requires scope `catalog:retrieve` or `catalog:retrieve-with-sequences`.

        ### Example Usage

        ```bash
        curl --request GET \
          --header "Authorization: Token ${KEY}" \
          --url https://api.developers.addgene.org/catalog/plasmid/
        ```
      parameters:
      - in: query
        name: article_authors
        schema:
          type: string
          minLength: 1
      - in: query
        name: article_pmid
        schema:
          type: string
          minLength: 1
      - in: query
        name: article_published
        schema:
          type: boolean
          nullable: true
      - in: query
        name: article_so
        schema:
          type: string
          minLength: 1
      - in: query
        name: article_title
        schema:
          type: string
          minLength: 1
      - in: query
        name: backbone
        schema:
          type: string
          minLength: 1
      - in: query
        name: bacterial_resistance
        schema:
          type: string
          minLength: 1
      - in: query
        name: catalog_item_id
        schema:
          type: string
          minLength: 1
      - in: query
        name: cloning_method
        schema:
          type: string
          minLength: 1
      - in: query
        name: experimental_use
        schema:
          type: string
          minLength: 1
      - in: query
        name: expression
        schema:
          type: string
          minLength: 1
      - in: query
        name: first_available_time
        schema:
          type: string
          format: date-time
      - in: query
        name: gene_ids
        schema:
          type: string
          minLength: 1
      - in: query
        name: genes
        schema:
          type: string
          minLength: 1
      - in: query
        name: is_industry
        schema:
          enum:
          - Available to Industry
          - null
        description: |-
          * `Available to Industry` - Available to Industry
          * `None` - None
      - in: query
        name: material_code
        schema:
          type: string
          minLength: 1
      - in: query
        name: mutations
        schema:
          type: string
          minLength: 1
      - in: query
        name: name
        schema:
          type: string
          minLength: 1
      - name: page
        required: false
        in: query
        description: A page number within the paginated result set.
        schema:
          type: integer
      - name: page_size
        required: false
        in: query
        description: Number of results to return per page.
        schema:
          type: integer
      - in: query
        name: pi_id
        schema:
          type: string
          minLength: 1
      - in: query
        name: pis
        schema:
          type: string
          minLength: 1
      - in: query
        name: plasmid_type
        schema:
          type: string
          minLength: 1
      - in: query
        name: promoters
        schema:
          type: string
          minLength: 1
      - in: query
        name: purpose
        schema:
          type: string
          minLength: 1
      - in: query
        name: resistance_marker
        schema:
          type: string
          minLength: 1
      - in: query
        name: sort_by
        schema:
          enum:
          - newest
          - id
          type: string
          default: id
          minLength: 1
        description: |-
          Sort by 'newest' or 'id'.

          * `newest` - Newest
          * `id` - ID
      - in: query
        name: species
        schema:
          type: string
          minLength: 1
      - in: query
        name: tags
        schema:
          type: string
          minLength: 1
      - in: query
        name: vector_types
        schema:
          type: string
          minLength: 1
      tags:
      - catalog
      security:
      - tokenAuth: []
      responses:
        '200':
          content:
            application/json:
              schema:
                $ref: '#/components/schemas/PaginatedPlasmidDocumentListList'
          description: ''

/catalog/plasmid/{id}/:
get:
operationId: catalog_plasmid_retrieve
description: |2

        ### Summary

        Endpoint to retrieve a single plasmid.

        Requires scope `catalog:retrieve` or `catalog:retrieve-with-sequences`.

        ### Example Usage

        ```bash
        curl --request GET \
          --header "Authorization: Token ${KEY}" \
          --url https://api.developers.addgene.org/catalog/plasmid/10878/
        ```
      parameters:
      - in: path
        name: id
        schema:
          type: integer
        description: Plasmid ID
        required: true
      tags:
      - catalog
      security:
      - tokenAuth: []
      responses:
        '200':
          content:
            application/json:
              schema:
                $ref: '#/components/schemas/Plasmid'
          description: ''

/catalog/plasmid-with-sequences/{id}/:
get:
operationId: catalog_plasmid_with_sequences_retrieve
description: |2

        ### Summary

        Endpoint to retrieve a single plasmid with associated sequences.

        Requires scope `catalog:retrieve-with-sequences`.

        ### Example Usage

        ```bash
        curl --request GET \
          --header "Authorization: Token ${KEY}" \
          --url https://api.developers.addgene.org/catalog/plasmid-with-sequences/10878/
        ```
      parameters:
      - in: path
        name: id
        schema:
          type: integer
        description: Plasmid ID
        required: true
      tags:
      - catalog
      security:
      - tokenAuth: []
      responses:
        '200':
          content:
            application/json:
              schema:
                $ref: '#/components/schemas/PlasmidWithSequences'
          description: ''

components:
schemas:
Article:
type: object
properties:
doi:
type: string
nullable: true
id:
type: integer
pubmed_id:
type: integer
nullable: true
url:
type: string
nullable: true
required: - doi - id - pubmed_id - url
Cloning:
type: object
properties:
backbone:
type: string
backbone_mutation:
type: string
backbone_origin:
type: string
backbone_size:
type: integer
nullable: true
promoter:
type: string
sequencing_primer_3:
type: string
sequencing_primer_5:
type: string
vector_types:
type: array
items:
type: string
required: - backbone - backbone_mutation - backbone_origin - backbone_size - promoter - sequencing_primer_3 - sequencing_primer_5
EntrezGene:
type: object
properties:
id:
type: integer
gene:
type: string
aliases:
type: string
required: - aliases - gene - id
Insert:
type: object
properties:
name:
type: string
alt_names:
type: array
items:
type: string
size:
type: integer
nullable: true
shRNA_sequence:
type: string
nullable: true
species:
type: array
items:
$ref: '#/components/schemas/Species'
genbank_ids:
type: array
items:
type: string
entrez_gene:
type: array
items:
$ref: '#/components/schemas/EntrezGene'
tags:
type: array
items:
$ref: '#/components/schemas/Tag'
mutation:
type: string
nullable: true
cloning:
allOf: - $ref: '#/components/schemas/InsertCloning'
nullable: true
required: - cloning - mutation - name - shRNA_sequence - size - species
InsertCloning:
type: object
properties:
clone_method:
type: string
promoter:
type: string
nullable: true
cloning_site_3:
type: string
nullable: true
cloning_site_5:
type: string
nullable: true
site_3_destroyed:
type: boolean
nullable: true
site_5_destroyed:
type: boolean
nullable: true
sequencing_primer_3:
type: string
nullable: true
sequencing_primer_5:
type: string
nullable: true
required: - clone_method - cloning_site_3 - cloning_site_5 - promoter - sequencing_primer_3 - sequencing_primer_5 - site_3_destroyed - site_5_destroyed
PaginatedPlasmidDocumentListList:
type: object
required: - count - results
properties:
count:
type: integer
example: 123
next:
type: string
nullable: true
format: uri
example: http://api.example.org/accounts/?page=4
previous:
type: string
nullable: true
format: uri
example: http://api.example.org/accounts/?page=2
results:
type: array
items:
$ref: '#/components/schemas/PlasmidDocumentList'
Plasmid:
type: object
properties:
id:
type: integer
viral_ids:
type: array
items:
type: string
name:
type: string
url:
type: string
description:
type: string
pi:
type: array
items:
type: integer
inserts:
type: array
items:
$ref: '#/components/schemas/Insert'
cloning:
$ref: '#/components/schemas/Cloning'
tags:
type: array
items:
$ref: '#/components/schemas/Tag'
bacterial_resistance:
type: string
growth_strain:
type: string
growth_temp:
type: integer
growth_notes:
type: string
plasmid_copy:
type: string
resistance_markers:
type: array
items:
type: string
origin:
type: string
depositor_comments:
type: string
article:
$ref: '#/components/schemas/Article'
terms:
type: string
required: - article - bacterial_resistance - cloning - depositor_comments - description - growth_notes - growth_strain - growth_temp - id - inserts - name - origin - plasmid_copy - tags - terms - url
PlasmidDocumentList:
type: object
description: Shared helpers for catalog document list serializers.
properties:
id:
type: integer
name:
type: string
purpose:
type: string
details_url:
type: string
format: uri
depositor:
type: array
items:
type: string
article:
type: array
items:
type: string
genes:
type: array
items:
type: string
species:
type: array
items:
type: string
experimental_use:
type: string
nullable: true
vector_types:
type: array
items:
type: string
expression:
type: array
items:
type: string
promoters:
type: array
items:
type: string
mutations:
type: array
items:
type: string
available_since:
type: string
format: date-time
nullable: true
required: - details_url - experimental_use - id - name - purpose
PlasmidWithSequences:
type: object
properties:
id:
type: integer
viral_ids:
type: array
items:
type: string
name:
type: string
url:
type: string
description:
type: string
pi:
type: array
items:
type: integer
inserts:
type: array
items:
$ref: '#/components/schemas/Insert'
cloning:
$ref: '#/components/schemas/Cloning'
tags:
type: array
items:
$ref: '#/components/schemas/Tag'
bacterial_resistance:
type: string
growth_strain:
type: string
growth_temp:
type: integer
growth_notes:
type: string
plasmid_copy:
type: string
resistance_markers:
type: array
items:
type: string
origin:
type: string
depositor_comments:
type: string
article:
$ref: '#/components/schemas/Article'
terms:
type: string
sequences:
$ref: '#/components/schemas/SequencesByType'
required: - article - bacterial_resistance - cloning - depositor_comments - description - growth_notes - growth_strain - growth_temp - id - inserts - name - origin - plasmid_copy - sequences - tags - terms - url
PlasmidsList:
type: object
properties:
plasmids:
type: array
items:
$ref: '#/components/schemas/Plasmid'
required: - plasmids
PlasmidsWithSequencesList:
type: object
properties:
plasmids:
type: array
items:
$ref: '#/components/schemas/PlasmidWithSequences'
required: - plasmids
Sequence:
type: object
properties:
sequence_description:
type: string
sequence:
type: string
genbank_url:
type: string
required: - genbank_url - sequence - sequence_description
SequencesByType:
type: object
properties:
public_user_full_sequences:
type: array
items:
$ref: '#/components/schemas/Sequence'
public_addgene_full_sequences:
type: array
items:
$ref: '#/components/schemas/Sequence'
public_user_partial_sequences:
type: array
items:
$ref: '#/components/schemas/Sequence'
public_addgene_partial_sequences:
type: array
items:
$ref: '#/components/schemas/Sequence'
required: - public_addgene_full_sequences - public_addgene_partial_sequences - public_user_full_sequences - public_user_partial_sequences
Species:
type: object
properties:
ncbi_taxonomy_id:
type: integer
nullable: true
species:
type: string
required: - ncbi_taxonomy_id - species
Tag:
type: object
properties:
tag:
type: string
location:
type: string
required: - location - tag
TestAuthResponse:
type: object
properties:
username:
type: string
email:
type: string
format: email
auth_method:
type: string
required: - auth_method - email - username
securitySchemes:
tokenAuth:
type: apiKey
in: header
name: Authorization
description: Token-based authentication with required prefix "Token"
servers:

- url: https://api.developers.addgene.org
  description: Addgene Developers API
  tags:
- name: auth
  description: |2

  ### Summary

  Authentication endpoints.

- name: download
  description: |2

  ### Summary

  Each download endpoint returns a 302 to a pre-signed S3 URL where the file can be downloaded directly.

  You will notice that the curl examples use the `--location` argument to follow the 302 redirects,
  which is necessary to reach the pre-signed S3 URL and ultimately receive a 200 for the download.
  If you are using a different HTTP client, please make sure it follows redirects as well.

  For security reasons, the pre-signed S3 URLs are short-lived. The default expiration time is 1 hour,
  which is subject to change without notice.
  If you receive an error when trying to download a file, please ensure that your api user has the correct
  scope and then try again to get a newly generated pre-signed URL.

- name: catalog
  description: |2

  ### Summary

  Catalog item endpoints.
