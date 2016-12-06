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

2. Set up the database tables with the schema in `db/schema.sql`:

   ````
   psql -h harps.guru -d orchestra -U $USERNAME

   orchestra=> \i db/schema.sql
   ````

3. Create or update `data/HARPS_all.csv`, which are the targets we will search 
   the ESO archive for.

4. From this working directory, run `scripts/eso_search_phase3.py`. This will 
   search ESO for Phase 3 products for every source (by `RA`, `Dec` position) in
   `data/HARPS_all.csv`:

   ````
   python scripts/eso_search_phase3.py
   ````

5. From this working directory, run `scripts/eso_retrieve.py`. This will request
   and download all data products (reduced spectra and intermediate data 
   products) identified in (4):
 
   ````
   python scripts/eso_retrieve.py
   ````

6. Follow the instructions printed at the end of the `scripts/eso_retrieve.py`
   script, namely:

   ````
   cd $DATA_DIR
   sh download_spectra.sh

   ls -lh *.tar | awk '{print "tar -xvf "" $9 "" --keep-old-files --force-local "}' > untar.sh
   sh untar.sh
   ````

7. From this working directory, run `db/ingest.py` to scrape the header 
   information from all spectra downloaded in (6), and ingest those headers into
   the database created in (1):

   ````
   python db/ingest.py
   ````


**Note**

If you need to run (5) multiple times, then you should run (6)-(7) before running
(5) again, because (5) will search for data products that have not been ingested.
