#!/usr/bin/env python3
"""Simple HTTP server for the pegamoid dissociation dashboard."""
import http.server
import os
import sys

PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 8080
WEB_DIR = os.path.dirname(os.path.abspath(__file__))


class Handler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=WEB_DIR, **kwargs)

    def log_message(self, fmt, *args):
        print(f"  {self.address_string()}  {fmt % args}")


print(f"\n  Pegamoid Dissociation Dashboard")
print(f"  Open in browser:  http://localhost:{PORT}")
print(f"  Serving from:     {WEB_DIR}")
print(f"  Press Ctrl+C to stop.\n")
with http.server.HTTPServer(("localhost", PORT), Handler) as httpd:
    httpd.serve_forever()
