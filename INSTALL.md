Setup
=====

1. Create a PostgreSQL database and put connection details (e.g., `hostname`, 
   `database`, `user`, `password`) in `db/credentials.yaml`. The contents
   of `db/credentials.yaml` should look something like:
   
   ````
   host: harps.guru
   database: orchestra
   user: $USERNAME
   password: $PASSWORD
   ````

2. Install [q3c](https://github.com/segasai/q3c) on the PostgreSQL database.


3. Set up the database tables with the schema in `db/schema.sql`:

   ````
   psql -h harps.guru -d orchestra -U $USERNAME

   orchestra=> \i db/schema.sql
   ````

4. Create or update `data/HARPS_all.csv`, which are the targets we will search 
   the ESO archive for.

5. From this working directory, run `scripts/eso_search_phase3.py`. This will 
   search ESO for Phase 3 products for every source (by `RA`, `Dec` position) in
   `data/HARPS_all.csv`:

   ````
   python scripts/eso_search_phase3.py
   ````

6. From this working directory, run `scripts/eso_retrieve.py`. This will request
   and download all data products (reduced spectra and intermediate data 
   products) identified in (5):
 
   ````
   python scripts/eso_retrieve.py
   ````

7. Follow the instructions printed at the end of the `scripts/eso_retrieve.py`
   script, namely:

   ````
   cd $DATA_DIR
   sh download_spectra.sh

   grep tar download_spectra.sh | awk '{print "tar -xvf " $1 " --keep-old-files --force-local "}' > untar.sh
   sh untar.sh
   ````

8. From this working directory, run `scripts/db_ingest_headers.py` to scrape the 
   header information from all spectra downloaded in (7), and ingest those headers
   into the database created in (1):

   ````
   python scripts/db_ingest_headers.py
   ````
 

**Note**

If you need to run (6) multiple times, then you should run (7)-(8) before running
(6) again, because (6) will search for data products that have not been ingested.
