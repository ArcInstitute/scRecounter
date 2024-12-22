# import
## batteries
import os
import warnings
from typing import List
## 3rd party
import psycopg2
import pandas as pd
from psycopg2.extras import execute_values
from psycopg2.extensions import connection
from tempfile import NamedTemporaryFile

# Suppress the specific warning
warnings.filterwarnings("ignore", message="pandas only supports SQLAlchemy connectable")

# functions
def db_connect() -> connection:
    """
    Connect to the sql database using SSL certificates.
    """
    # get certs
    certs = get_db_certs()
    database = "sragent-test" if os.environ["GCP_SQL_DB_TENANT"].upper() == "TEST" else "sragent"
    # connect
    db_params = {
        'host': "35.243.133.29", #os.environ["GCP_SQL_DB_HOST"],
        'database': database, # os.environ["GCP_SQL_DB_NAME"],
        'user': "postgres", #os.environ["GCP_SQL_DB_USERNAME"],
        'password': os.getenv("GCP_SQL_DB_PASSWORD", get_secret("GCP_SQL_DB_PASSWORD")),
        'sslmode': 'verify-ca',
        'sslrootcert': certs["server-ca.pem"],
        'sslcert': certs["client-cert.pem"],
        'sslkey': certs["client-key.pem"],
        'port': '5432',
        'connect_timeout': 30
    }
    conn = psycopg2.connect(**db_params)
    # delete certs
    for cert in certs.values():
        os.remove(cert)
    return conn

def add_to_log(
        df, sample: str, accession: str, process: str, step: str, status: str, msg: str
        ) -> pd.DataFrame:
    """
    Add log entry to dataframe.
    Args:
        log_df: Log dataframe
        sample: Sample name
        accession: SRA accession
        process: Process name
        step: Step name
        status: Status
        msg: Message
    Returns:
        pd.DataFrame: Updated log dataframe
    """
    if len(msg) > 200:
        msg = str(msg[:(200-3)]) + '...'
    df.loc[len(df)] = [sample, accession, process, step, status, msg]

def db_upsert(df: pd.DataFrame, table_name: str, conn: connection) -> None:
    """
    Upload a pandas DataFrame to PostgreSQL, performing an upsert operation.
    If records exist (based on unique constraints), update them; otherwise insert new records.
    Args:
        df: pandas DataFrame to upload
        table_name: name of the target table
        conn: psycopg2 connection object
    """   
    # if df is empty, return
    if df.empty:
        return
    # if df is not dataframe, try to convert
    if not isinstance(df, pd.DataFrame):
        try:
            df = pd.DataFrame(df)
        except Exception as e:
            raise Exception(f"Error converting input to DataFrame: {str(e)}")

    # filter to overlapping target columns
    df = df[list(set(get_table_columns(table_name, conn)).intersection(df.columns))]

    # Get DataFrame columns
    columns = list(df.columns)
    
    # Create ON CONFLICT clause based on unique constraints
    unique_columns = get_unique_columns(table_name, conn)

    # Exclude 'id' column from the upsert
    if "id" in columns:
        df = df.drop(columns=["id"])
        columns.remove("id")

    # Drop duplicate records based on unique columns
    df.drop_duplicates(subset=unique_columns, keep='first').copy()

    # Convert DataFrame to list of tuples
    values = [tuple(x) for x in df.to_numpy()]

    # Create the INSERT statement with ON CONFLICT clause
    insert_stmt = f"INSERT INTO {table_name} ({', '.join(columns)})"
    insert_stmt += f"\nVALUES %s"

    # Add DO UPDATE SET clause for non-unique columns
    do_update_set = [col for col in columns if col not in unique_columns]
    if do_update_set:
        do_update_set = ', '.join(f"{col} = EXCLUDED.{col}" for col in do_update_set)
        insert_stmt += f"\nON CONFLICT ({', '.join(unique_columns)})"
        insert_stmt += f"\nDO UPDATE SET {do_update_set}"
    else:
        # if no non-unique columns, add DO NOTHING clause
        insert_stmt += f"\nON CONFLICT ({', '.join(unique_columns)}) DO NOTHING"

    # Execute the query
    try:
        with conn.cursor() as cur:
            execute_values(cur, insert_stmt, values)
            conn.commit()
    except Exception as e:
        conn.rollback()
        raise Exception(f"Error uploading data to {table_name}: {str(e)}")

def get_table_columns(table: str, conn: connection) -> List[str]:
    """
    Get column names for a table from the database schema.
    Args:
        table: Name of the table
        conn: Database connection 
    Returns:
        List of column names
    """
    query = """
    SELECT column_name
    FROM information_schema.columns
    WHERE table_name = %s;
    """
    
    with conn.cursor() as cur:
        cur.execute(query, (table,))
        columns = cur.fetchall()
    return [col[0] for col in columns]

def get_unique_columns(table: str, conn: connection) -> List[str]:
    """
    Get all unique constraint columns for a table from the database schema.
    Prioritizes composite unique constraints over primary keys.
    Args:
        table: Name of the table
        conn: Database connection 
    Returns:
        List of column names that form the most appropriate unique constraint
    """
    query = """
    SELECT c.contype, ARRAY_AGG(a.attname ORDER BY array_position(c.conkey, a.attnum)) as columns
    FROM pg_constraint c
    JOIN pg_class t ON c.conrelid = t.oid
    JOIN pg_attribute a ON a.attrelid = t.oid AND a.attnum = ANY(c.conkey)
    WHERE t.relname = %s 
    AND c.contype IN ('p', 'u')  -- primary key or unique constraint
    GROUP BY c.conname, c.contype
    ORDER BY c.contype DESC;  -- 'u'nique before 'p'rimary key
    """
    
    with conn.cursor() as cur:
        cur.execute(query, (table,))
        constraints = cur.fetchall()
        
    if not constraints:
        raise ValueError(f"No unique constraints found in table {table}")
    
    # Prefer composite unique constraints over single-column primary keys
    for constraint_type, columns in constraints:
        if len(columns) > 1 or constraint_type == 'u':
            return columns
    
    # Fall back to primary key if no other suitable constraint found
    return constraints[0][1]

def get_secret(secret_id: str) -> str:
    """
    Fetch secret from GCP Secret Manager.
    Rquired environment variables: GCP_PROJECT_ID, GOOGLE_APPLICATION_CREDENTIALS
    Args:
        secret_id: The secret id
    Returns:
        The secret value
    """
    from google.auth import default
    from google.cloud import secretmanager

    _, project_id = default()  # Use default credentials; project_id is inferred
    if not project_id:
        project_id = os.environ["GCP_PROJECT_ID"]
    name = f"projects/{project_id}/secrets/{secret_id}/versions/latest"
    client = secretmanager.SecretManagerServiceClient()
    response = client.access_secret_version(request={"name": name})
    return response.payload.data.decode('UTF-8')

def get_db_certs(certs=["server-ca.pem", "client-cert.pem", "client-key.pem"]) -> dict:
    """
    Download certificates from GCP Secret Manager and save them to temporary files.
    Args:
        certs: A list of certificate ids
    Returns:
        A dictionary containing the paths to the temporary files
    """
    idx = {
        "server-ca.pem": "SRAgent_db_server_ca",
        "client-cert.pem": "SRAgent_db_client_cert",
        "client-key.pem": "SRAgent_db_client_key"
    }
    cert_files = {}
    for cert in certs:
        cert_files[cert] = download_secret(idx[cert])
    return cert_files

def download_secret(secret_id: str) -> str:
    """
    Download a secret from GCP Secret Manager and save it to a temporary file.
    Args:
        secret_id: The secret id
    Returns:
        The path to the temporary file containing the secret
    """
    secret_value = get_secret(secret_id)
    temp_file = NamedTemporaryFile(delete=False, mode='w', encoding='utf-8')
    with temp_file as f:
        f.write(secret_value)
        f.flush()
    return temp_file.name

def get_srx_metadata_limit5(conn):
    query = """
    SELECT * FROM srx_metadata LIMIT 5;
    """
    return pd.read_sql(query, conn)

# main
if __name__ == "__main__":
    from dotenv import load_dotenv
    load_dotenv()
    #print(get_db_certs())
    with db_connect() as conn:
        print(get_srx_metadata_limit5(conn))