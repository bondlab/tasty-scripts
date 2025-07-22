#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 07:26:03 2025
GUI for running FeatureBot script
"""

import sys
import os
from PySide6.QtWidgets import (
    QApplication, QWidget, QLabel, QLineEdit, QPushButton, QFileDialog,
    QVBoxLayout, QHBoxLayout, QCheckBox, QGroupBox, QRadioButton, QTextEdit
)
from PySide6.QtCore import Qt
import subprocess

class FeatureBotGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("FeatureBot")
        self.setMinimumWidth(500)

        layout = QVBoxLayout()

        # GenBank file or directory
        self.gbk_path = QLineEdit()
        self.gbk_path.setToolTip("Select a single GenBank file or a folder for batch mode.")
        gbk_path_layout = QHBoxLayout()
        gbk_path_layout.addWidget(self.gbk_path)

        gbk_file_btn = QPushButton("📄 File")
        gbk_file_btn.clicked.connect(self.select_gbk_file)
        gbk_path_layout.addWidget(gbk_file_btn)

        gbk_folder_btn = QPushButton("📁 Folder")
        gbk_folder_btn.clicked.connect(self.select_gbk_folder)
        gbk_path_layout.addWidget(gbk_folder_btn)

        layout.addWidget(QLabel("GenBank input:"))
        layout.addLayout(gbk_path_layout)

        # Feature file or directory
        self.feature_path = QLineEdit()
        self.feature_path.setToolTip("Select a feature file (GFF, TSV) or a folder for batch mode.")
        feat_path_layout = QHBoxLayout()
        feat_path_layout.addWidget(self.feature_path)

        feat_file_btn = QPushButton("📄 File")
        feat_file_btn.clicked.connect(self.select_feat_file)
        feat_path_layout.addWidget(feat_file_btn)

        feat_folder_btn = QPushButton("📁 Folder")
        feat_folder_btn.clicked.connect(self.select_feat_folder)
        feat_path_layout.addWidget(feat_folder_btn)

        layout.addWidget(QLabel("Feature file or folder:"))
        layout.addLayout(feat_path_layout)

        # Mode selection
        mode_group = QGroupBox("Annotation mode")
        mode_layout = QVBoxLayout()
        self.signal_btn = QRadioButton("Signal peptide analysis (SignalP 6.0)")
        self.signal_btn.setToolTip("Add SignalP peptide annotations from a GFF file.")
        self.beta_btn = QRadioButton("Beta barrel analysis (TMHMM)")
        self.multiheme_btn = QRadioButton("Multiheme cytochromes (hemehunter)")
        self.cluster_btn = QRadioButton("Clusters (wirefinder, etc)")
        mode_layout.addWidget(self.signal_btn)
        mode_layout.addWidget(self.beta_btn)
        mode_layout.addWidget(self.multiheme_btn)
        mode_layout.addWidget(self.cluster_btn)
        mode_group.setLayout(mode_layout)
        layout.addWidget(mode_group)

        # Optional flags
        self.force_box = QCheckBox("Overwrite prior analyses (--force)")
        self.force_box.setToolTip("If checked, existing output files will be overwritten.")
        self.proteinid_box = QCheckBox("Match protein ID instead of locus tag (--proteinid)")
        layout.addWidget(self.force_box)
        layout.addWidget(self.proteinid_box)

        # Run button
        run_button = QPushButton("Run FeatureBot")
        run_button.setMinimumHeight(40)
        run_button.clicked.connect(self.run_featurebot)
        layout.addWidget(run_button)

        # Output area
        self.output = QTextEdit()
        self.output.setReadOnly(True)
        layout.addWidget(self.output)

        self.setLayout(layout)

    def select_gbk_file(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select GenBank file", "", "GenBank Files (*.gb *.gbk *.gbff)")
        if file:
            self.gbk_path.setText(file)

    def select_gbk_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Select GenBank folder")
        if folder:
            self.gbk_path.setText(folder)

    def select_feat_file(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select feature file", "", "All Files (*.*)")
        if file:
            self.feature_path.setText(file)

    def select_feat_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Select feature folder")
        if folder:
            self.feature_path.setText(folder)

    def run_featurebot(self):
        script_path = os.path.join(os.path.dirname(__file__), "featurebot.py")
        if not os.path.isfile(script_path):
            self.output.append(f"❌ featurebot.py not found at {script_path}")
            return
        args = [sys.executable, script_path]

        gbk_input = self.gbk_path.text()
        feat_input = self.feature_path.text()

        if not os.path.exists(gbk_input):
            self.output.append(f"❌ GenBank path does not exist: {gbk_input}")
            return

        if not os.path.exists(feat_input):
            self.output.append(f"❌ Feature path does not exist: {feat_input}")
            return

        if os.path.isdir(gbk_input):
            args += ["-gbkdir", gbk_input]
        else:
            args += ["-gbk", gbk_input]

        if os.path.isdir(feat_input):
            args += ["-featuredir", feat_input]
        else:
            args += ["-gff", feat_input]

        if self.signal_btn.isChecked():
            args.append("-signal")
        elif self.beta_btn.isChecked():
            args.append("-beta")
        elif self.multiheme_btn.isChecked():
            args.append("-multiheme")
        elif self.cluster_btn.isChecked():
            args.append("-cluster")
        else:
            self.output.append("❌ Please select one annotation mode.")
            return

        if self.force_box.isChecked():
            args.append("--force")
        if self.proteinid_box.isChecked():
            args.append("--proteinid")

        self.output.append(f"➡️ Running: {' '.join(args)}\n")

        try:
            result = subprocess.run(args, text=True, capture_output=True)
            self.output.append(result.stdout)
            if result.stderr:
                self.output.append("\n[stderr]\n" + result.stderr)
        except Exception as e:
            self.output.append(f"❌ Error: {str(e)}")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = FeatureBotGUI()
    window.show()
    sys.exit(app.exec())
