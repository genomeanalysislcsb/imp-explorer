# Run the app

```bash
imp-explorer run --data release_v0.1                 # http://127.0.0.1:8050
imp-explorer run --data release_v0.1 --port 8765     # another port
```

| option | default | meaning |
|---|---|---|
| `--data` | – | release directory |
| `--host` | `127.0.0.1` | `0.0.0.0` to accept connections from other machines (needed for an SSH tunnel to a compute node) |
| `--port` | 8050 | port |
| `--plugin` | – | module with extra pages (repeatable), see [plug-ins](plugins.md) |
| `--table` | – | extra release table a plug-in needs (repeatable) |
| `--debug` | off | Dash debug mode (auto-reload, error pop-ups) |

## On an HPC compute node

Start the app on the node, then forward the port through the login node from your laptop:

```bash
# on the compute node (e.g. aion-0115)
imp-explorer run --data release_v0.1 --host 0.0.0.0 --port 8050

# on your laptop (ULHPC login nodes listen on port 8022)
ssh -p 8022 -L 8050:aion-0115:8050 <user>@access-aion.uni.lu
```

Open <http://localhost:8050>. The forwarding is `local_port:node:remote_port`; all three parts are needed.

## On a laptop

Copy the release directory (Parquet files + manifest) and run the app there; nothing else is needed.

## In production

`imp-explorer run` uses Flask's development server, which is fine for a few users. For a public deployment run the
Dash app behind a WSGI server:

```python
# wsgi.py
from impexplorer.app import create_app
app = create_app("/data/release_v1.0")
server = app.server
```

```bash
gunicorn -w 4 -b 0.0.0.0:8050 wsgi:server
```
