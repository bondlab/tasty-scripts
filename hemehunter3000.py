#!/usr/bin/env python3

import sys
import os
import io
import contextlib
from PyQt5.QtWidgets import (
    QApplication, QWidget, QLabel, QPushButton,
    QVBoxLayout, QHBoxLayout, QFileDialog, QLineEdit,
    QCheckBox, QSpinBox, QTextEdit
)
from PyQt5.QtGui import QTextCursor

# Import your script's functions directly
import hemehunter

# === Helper function to get help text ===
def get_help_text():
    parser = hemehunter.build_parser()  # This makes the help work
    help_buf = io.StringIO()
    with contextlib.redirect_stdout(help_buf):
        parser.print_help()
    return help_buf.getvalue()


class HemeHunterGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("HemeHunter 3000")

        # === Input file or folder ===
        self.input_label = QLabel("GenBank file or folder: TSV output will be generated for each")
        self.input_path = QLineEdit()

        self.browse_button = QPushButton("Browse files")
        self.browse_button.clicked.connect(self.browse_file)

        self.browse_folder_button = QPushButton("Browse folders")
        self.browse_folder_button.clicked.connect(self.browse_folder)

        # Group input path + browse buttons in a horizontal row
        path_layout = QHBoxLayout()
        path_layout.addWidget(self.input_path)
        path_layout.addWidget(self.browse_button)
        path_layout.addWidget(self.browse_folder_button)

        # === Options ===
        self.fasta_checkbox = QCheckBox("Also output FASTA of all cytochromes")
       
        self.force_checkbox = QCheckBox("Force overwrite of previous output")
        
        self.translations_checkbox = QCheckBox("Translate missing CDS if missing from GenBank file")

        self.motif_label = QLabel("Mininmim number of motifs:")
        self.motif_spin = QSpinBox()
        self.motif_spin.setMinimum(1)
        self.motif_spin.setValue(3)

        # === Help button ===
        self.help_button = QPushButton("Show Help")
        self.help_button.clicked.connect(self.show_help)

        # === Run button ===
        self.run_button = QPushButton("-- Hunt for some Hemes --")
        self.run_button.clicked.connect(self.run_hemehunter)
        # Make it bigger
        self.run_button.setFixedHeight(50)  # height
        font = self.run_button.font()
        font.setPointSize(14)               # bigger text
        self.run_button.setFont(font)
        
        # === Log output ===
        self.log_output = QTextEdit()
        self.log_output.setReadOnly(True)

        # === Main layout ===
        layout = QVBoxLayout()
        layout.addWidget(self.input_label)
        layout.addLayout(path_layout)
        layout.addWidget(self.fasta_checkbox)
        layout.addWidget(self.translations_checkbox)
        layout.addWidget(self.force_checkbox)
        layout.addWidget(self.motif_label)
        layout.addWidget(self.motif_spin)
        layout.addWidget(self.help_button)
        layout.addWidget(self.run_button)
        layout.addWidget(QLabel("Log Output:"))
        layout.addWidget(self.log_output)

        self.setLayout(layout)

        # Optional: set initial window size
        self.resize(800, 600)

    def show_help(self):
        help_text = get_help_text()
        self.log_output.clear()
        self.log_output.append(help_text)

    def browse_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select GenBank file",
            "",
            "GenBank Files (*.gb *.gbk *.gbff);;All Files (*)"
        )
        if path:
            self.input_path.setText(path)

    def browse_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Select a Folder")
        if folder:
            self.input_path.setText(folder)

    def run_hemehunter(self):
        # Mimic argparse.Namespace with a simple Args class
        class Args:
            genbank_file = self.input_path.text()
            output = None  # Output file path not used in this minimal GUI
            fasta = self.fasta_checkbox.isChecked()
            translations = self.translations_checkbox.isChecked()
            force = self.force_checkbox.isChecked()
            min_motifs = self.motif_spin.value()

        if not Args.genbank_file or not os.path.exists(Args.genbank_file):
            self.log_output.append("[ERROR] Invalid input file or folder.")
            return

        # Redirect stdout to the GUI log box
        self.log_output.clear()
        old_stdout = sys.stdout
        sys.stdout = self

        try:
            if os.path.isdir(Args.genbank_file):
                self.log_output.append(f"[INFO] Batch mode for folder: {Args.genbank_file}")
                for fname in os.listdir(Args.genbank_file):
                    if fname.endswith((".gb", ".gbk", ".gbff")):
                        full_path = os.path.join(Args.genbank_file, fname)
                        print(f"📄 Processing: {fname}")
                        hemehunter.run_on_one_genbank_file(full_path, Args)
            else:
                hemehunter.run_on_one_genbank_file(Args.genbank_file, Args)
            self.log_output.append("\n[INFO] Finished processing.")
        except Exception as e:
            self.log_output.append(f"[ERROR] {e}")
        finally:
            sys.stdout = old_stdout

    # To make print() output appear in the log_output text box
    def write(self, message):
        self.log_output.moveCursor(QTextCursor.End)
        self.log_output.insertPlainText(message)
        self.log_output.moveCursor(QTextCursor.End)

    def flush(self):
        pass  # Required by Python’s print system


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = HemeHunterGUI()
    window.show()
    sys.exit(app.exec_())