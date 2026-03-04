import os
import json
from pathlib import Path

def main():
    ui_dir = Path(__file__).parent
    data_dir = ui_dir / "data"
    output_js = ui_dir / "data.js"

    results = {}
    
    if not data_dir.exists():
        print(f"Error: {data_dir} does not exist.")
        return

    print(f"Scanning {data_dir} for MAGICAL outputs...")
    
    for file_path in data_dir.glob("*.txt"):
        name = file_path.stem
        try:
            with open(file_path, 'r') as f:
                content = f.read()
                results[name] = content
            print(f"  Added: {name}")
        except Exception as e:
            print(f"  Failed to read {name}: {e}")

    if not results:
        print("No .txt files found in data folder.")
        return

    # Write to a JS file that defines a global object
    # Using JSON.stringify makes it safe for embedding
    js_content = f"window.MAGICAL_RESULTS = {json.dumps(results, indent=2)};"
    
    with open(output_js, 'w') as f:
        f.write(js_content)
    
    print(f"\nSuccessfully generated {output_js}")
    print("Open magical_explorer.html to view results.")

if __name__ == "__main__":
    main()
