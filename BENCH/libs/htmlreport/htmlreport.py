##########################################################
#  CROCO under CeCILL-C
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# python
import os
import json
import shutil

# internal
from ..helpers import Messaging


def format_status(value):
    if value is True:
        format_td = '<td class="status-true">Done</td>'
    elif value is False:
        format_td = '<td class="status-false">Error</td>'
    else:
        format_td = '<td class="status-null">N/A</td>'
    return format_td


def generate_status_html(data):
    html_rows = ""
    for key, entry in data.items():
        case = entry["case"]
        if entry["restarted"]:
            variant = entry["variant"] + "-rst"
        else:
            variant = entry["variant"]
        status = entry["status"]
        html_rows += f"""
            <tr>
               <td>{case}</td>
               <td>{variant}</td>
               {format_status(status["build"]["status"])}
               {format_status(status["run"]["status"])}
               {format_status(status["check"]["status"])}
               {format_status(status["plotphy"]["status"])}
            </tr>
            """
    html_table = f"""
        <table>
            <thead>
                  <tr>
                     <th>Case</th>
                     <th>Variant</th>
                     <th>Build Status</th>
                     <th>Run Status</th>
                     <th>Check Status</th>
                     <th>Plotphy Status</th>
                  </tr>
            </thead>
            <tbody>
                  {html_rows}
            </tbody>
        </table>"""

    return html_table


def generate_treeview_json(root_path):
    def traverse_directory(path, base_url=""):
        nodes = []
        files = []
        directories = []
        status_report_node = None  # To store the status_report.json node

        for name in sorted(os.listdir(path), key=str.lower):
            full_path = os.path.join(path, name)
            relative_path = os.path.relpath(full_path, root_path)
            url_path = os.path.join(base_url, relative_path).replace(os.sep, "/")

            if os.path.isdir(full_path):
                if not os.path.islink(full_path):
                    directories.append(
                        {"text": name, "nodes": traverse_directory(full_path, base_url)}
                    )
            else:
                if name == "status_report.json":
                    # Special handling for status_report.json
                    status_report_node = {
                        "text": generate_status_html(json.load(open(full_path, "r")))
                    }

                elif ("treeview.html" not in url_path) and (
                    "styles.css" not in url_path
                ):
                    files.append(
                        {"text": f'<a href="{url_path}" target="main-frame">{name}</a>'}
                    )

        # Sort files and directories alphabetically
        files = sorted(files, key=lambda x: x["text"].lower())
        directories = sorted(directories, key=lambda x: x["text"].lower())

        # Place status_report.json node first (if it exists), then sorted files, then sorted directories
        if status_report_node:
            nodes.append(status_report_node)
        nodes.extend(files)
        nodes.extend(directories)

        return nodes

    tree_structure = traverse_directory(root_path)

    return tree_structure


# Function to convert JSON to HTML
def json_to_html(data):
    html_content = ""

    for item in data:
        if isinstance(item, dict):
            text = item.get("text", "")
            # Create a collapsible section using <details> and <summary>
            if "nodes" in item:
                html_content += "<details>\n"
                html_content += f"  <summary>{text}</summary>\n"
                html_content += '  <div class="tree">\n'
                html_content += json_to_html(item["nodes"])  # Recursively process nodes
                html_content += "  </div>\n"
                html_content += "</details>\n"
            else:
                html_content += text + "\n"
    return html_content


# Function to convert JSON tree to HTML
def generate_html(
    data_path,
    output_file="treeview.html",
    navigation_content="",
    style_file="libs/htmlreport/styles.css",
):
    json_data = generate_treeview_json(data_path)

    html_output = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CROCO test report</title>
    <link rel="stylesheet" href="styles.css">
</head>
<body>
    %s
    <main>
    <h1>CROCO test report - %s</h1>

    <div class="tree">
""" % (
        navigation_content,
        os.path.split(data_path)[-1],
    )

    html_output += json_to_html(json_data)

    html_output += """
    </div>

    <iframe id="main-frame" name="main-frame"></iframe>
    </main>
</body>
</html>
"""

    # Write the generated HTML to a file
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(html_output)
    # copy css file
    shutil.copyfile(
        style_file,
        os.path.join(data_path, "styles.css"),
    )

    Messaging.step(f"HTML file generated: {output_file}")


# Function to create index html file for all results
def generate_global_html(base_dir, style_file="libs/htmlreport/styles.css"):
    output_file = os.path.join(base_dir, "index.html")

    # Find all folders in base_dir
    folders = [
        folder
        for folder in os.listdir(base_dir)
        if os.path.isdir(os.path.join(base_dir, folder))
        and not os.path.islink(os.path.join(base_dir, folder))
    ]

    # navigation in folders (relative paths using ../)
    nav_in_folder = f"""  <nav>
      <h2><a href="../index.html">Results</a></h2>
      <ul>
         {
        "".join(
            f'<li><a href="../{folder}/treeview.html">View {folder}</a></li>'
            for folder in folders
        )
    }
      </ul>
      </nav>"""

    for folder in folders:
        if not os.path.islink(os.path.join(base_dir, folder)):
            generate_html(
                os.path.join(base_dir, folder),
                output_file=os.path.join(
                    os.path.join(base_dir, folder), "treeview.html"
                ),
                navigation_content=nav_in_folder,
            )

    # HTML content
    html_content = f"""<!DOCTYPE html>
   <html lang="en">
   <head>
   <meta charset="UTF-8">
   <meta name="viewport" content="width=device-width, initial-scale=1.0">
   <title>CROCO test report</title>
   <link rel="stylesheet" href="styles.css">
   </head>
   <body>
   <nav>
      <h2><a href="index.html">Results</a></h2>
      <ul>
         {
        "".join(
            f'<li><a href="{folder}/treeview.html">View {folder}</a></li>'
            for folder in folders
        )
    }
      </ul>
   </nav>
   <main>
      <h1>CROCO test report</h1>
      <p>Use the navigation menu on the left to explore the treeviews in each folder.</p>
   </main>
   </body>
   </html>
   """

    # write HTML index
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(html_content)
    # copy css file
    shutil.copyfile(
        style_file,
        os.path.join(base_dir, "styles.css"),
    )

    Messaging.step(f"Index file generated: {os.path.abspath(output_file)}")
