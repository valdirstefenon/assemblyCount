import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import os
import gzip
import csv
import threading

# Global configurations for Matplotlib (Publication quality)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2


class GenomeAnalyzerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("assemblyCount")
        self.root.geometry("1050x850")
        self.root.configure(bg="#ecf0f1")

        # ttk Styling (Blue tones)
        self.style = ttk.Style()
        self.style.theme_use('clam')

        self.bg_color = "#ecf0f1"
        self.primary_blue = "#2980b9"
        self.dark_blue = "#1f618d"
        self.text_color = "#2c3e50"

        self.style.configure('TFrame', background=self.bg_color)
        self.style.configure('TLabel', background=self.bg_color, foreground=self.text_color, font=('Arial', 10))
        self.style.configure('TButton', background=self.primary_blue, foreground='white', font=('Arial', 10, 'bold'), borderwidth=0)
        self.style.map('TButton', background=[('active', self.dark_blue)])
        self.style.configure('Header.TLabel', font=('Arial', 14, 'bold'), foreground=self.dark_blue)
        self.style.configure('TLabelframe', background=self.bg_color, foreground=self.dark_blue, font=('Arial', 10, 'bold'))
        self.style.configure('TLabelframe.Label', background=self.bg_color, foreground=self.dark_blue)
        self.style.configure('TCheckbutton', background=self.bg_color, foreground=self.text_color)

        # State attributes
        self.filepath = None
        self.names = []
        self.lengths = []
        self.total_length = 0
        self.gc_content = 0.0
        self.n50, self.l50 = 0, 0
        self.n75, self.l75 = 0, 0
        self.n90, self.l90 = 0, 0
        self.aux = 0.0
        self.ambiguous_counts = {}
        self.min_length = 0

        self.create_widgets()

    def create_widgets(self):
        # Top Frame (Buttons and Status)
        top_frame = ttk.Frame(self.root, padding="10")
        top_frame.pack(fill=tk.X)

        ttk.Label(top_frame, text="assemblyCount", style='Header.TLabel').pack(side=tk.LEFT, padx=10)

        self.btn_load = ttk.Button(top_frame, text="Load FASTA", command=self.load_fasta)
        self.btn_load.pack(side=tk.RIGHT, padx=5)

        self.btn_process = ttk.Button(top_frame, text="Process & Save Results", command=self.process_data_thread, state=tk.DISABLED)
        self.btn_process.pack(side=tk.RIGHT, padx=5)

        # Options Frame (Minimum length filter + gz support note)
        options_frame = ttk.LabelFrame(self.root, text=" Options ", padding="10")
        options_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(options_frame, text="Minimum scaffold length (bp):").pack(side=tk.LEFT, padx=5)
        self.min_length_var = tk.StringVar(value="0")
        self.entry_min_length = ttk.Entry(options_frame, textvariable=self.min_length_var, width=12)
        self.entry_min_length.pack(side=tk.LEFT, padx=5)

        ttk.Label(options_frame, text="(0 = no filter; common: 1000)").pack(side=tk.LEFT, padx=5)

        # Middle Frame (Text Information)
        middle_frame = ttk.Frame(self.root, padding="10")
        middle_frame.pack(fill=tk.BOTH, expand=False)

        info_frame = ttk.LabelFrame(middle_frame, text=" Assembly Results ", padding="10")
        info_frame.pack(fill=tk.X)

        self.text_output = tk.Text(info_frame, height=18, width=80, bg="white", fg=self.text_color,
                                   font=('Courier', 10), relief="flat", borderwidth=1)
        self.text_output.pack(fill=tk.X)

        # Status bar
        self.status_var = tk.StringVar(value="Ready.")
        self.status_bar = ttk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.pack(fill=tk.X, side=tk.BOTTOM)

        # Bottom Frame (Plot)
        self.plot_frame = ttk.LabelFrame(self.root, text=" Cumulative Bar Chart ", padding="10")
        self.plot_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

    # ------------------------------------------------------------------
    # File loading
    # ------------------------------------------------------------------
    def load_fasta(self):
        self.filepath = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA Files", "*.fasta *.fa *.fna *.fas *.fasta.gz *.fa.gz *.fna.gz *.fas.gz"),
                       ("All Files", "*.*")]
        )
        if self.filepath:
            self.btn_process.config(state=tk.NORMAL)
            self.text_output.delete(1.0, tk.END)
            self.text_output.insert(tk.END, f"File loaded: {os.path.basename(self.filepath)}\n")
            self.text_output.insert(tk.END, "Click 'Process & Save Results' to start the analysis.\n")
            self.status_var.set(f"Loaded: {os.path.basename(self.filepath)}")

    # ------------------------------------------------------------------
    # FASTA parsing (supports .gz)
    # ------------------------------------------------------------------
    def parse_fasta(self):
        """Reads FASTA (optionally gzipped) and returns names, lengths, GC, AT and ambiguous counts."""
        names = []
        lengths = []
        gc_count = 0
        at_count = 0
        ambiguous_counts = {}

        current_len = 0
        current_name = None

        if self.filepath.endswith('.gz'):
            f = gzip.open(self.filepath, 'rt')
        else:
            f = open(self.filepath, 'r')

        with f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_name is not None:
                        names.append(current_name)
                        lengths.append(current_len)
                    current_name = line[1:].split()[0]
                    current_len = 0
                else:
                    seq = line.upper()
                    current_len += len(seq)
                    gc_count += seq.count('G') + seq.count('C')
                    at_count += seq.count('A') + seq.count('T')
                    for ch in seq:
                        if ch not in ('A', 'T', 'G', 'C'):
                            ambiguous_counts[ch] = ambiguous_counts.get(ch, 0) + 1
            if current_name is not None:
                names.append(current_name)
                lengths.append(current_len)

        # Sort by length descending, keeping name pairs
        paired = sorted(zip(names, lengths), key=lambda x: x[1], reverse=True)
        names_sorted = [p[0] for p in paired]
        lengths_sorted = [p[1] for p in paired]

        return names_sorted, lengths_sorted, gc_count, at_count, ambiguous_counts

    # ------------------------------------------------------------------
    # Nx / Lx calculation
    # ------------------------------------------------------------------
    def calculate_nx_lx(self, lengths, total_len, x):
        """Calculates Nx and Lx (e.g., N50 and L50)."""
        if total_len == 0:
            return 0, 0
        target = total_len * (x / 100.0)
        cum_sum = 0
        for i, length in enumerate(lengths):
            cum_sum += length
            if cum_sum >= target:
                return length, (i + 1)
        return 0, 0

    # ------------------------------------------------------------------
    # Threaded processing (avoids freezing the GUI)
    # ------------------------------------------------------------------
    def process_data_thread(self):
        self.btn_process.config(state=tk.DISABLED)
        self.btn_load.config(state=tk.DISABLED)
        self.status_var.set("Processing... please wait.")
        thread = threading.Thread(target=self.process_data)
        thread.daemon = True
        thread.start()

    def process_data(self):
        if not self.filepath:
            return

        try:
            # Validate min_length input
            try:
                self.min_length = int(self.min_length_var.get())
                if self.min_length < 0:
                    raise ValueError
            except ValueError:
                self.root.after(0, lambda: messagebox.showerror("Error", "Minimum length must be a non-negative integer."))
                self.reset_buttons()
                return

            self.root.after(0, lambda: self.text_output.delete(1.0, tk.END))
            self.root.after(0, lambda: self.text_output.insert(tk.END, "Processing data...\n"))
            self.root.after(0, lambda: self.status_var.set("Reading FASTA file..."))

            # 1. Parse FASTA
            names, lengths, gc_count, at_count, ambiguous_counts = self.parse_fasta()

            if not lengths:
                self.root.after(0, lambda: messagebox.showerror("Error", "No sequences found in the FASTA file."))
                self.reset_buttons()
                return

            total_before_filter = len(lengths)
            total_size_before_filter = sum(lengths)

            # 2. Apply minimum length filter
            filtered = [(n, l) for n, l in zip(names, lengths) if l >= self.min_length]
            if not filtered:
                self.root.after(0, lambda: messagebox.showerror("Error", "No scaffolds remain after applying the minimum length filter."))
                self.reset_buttons()
                return

            self.names = [x[0] for x in filtered]
            self.lengths = [x[1] for x in filtered]

            # 3. Total assembly size (after filter)
            self.total_length = sum(self.lengths)

            # 4. Nx / Lx
            self.n50, self.l50 = self.calculate_nx_lx(self.lengths, self.total_length, 50)
            self.n75, self.l75 = self.calculate_nx_lx(self.lengths, self.total_length, 75)
            self.n90, self.l90 = self.calculate_nx_lx(self.lengths, self.total_length, 90)

            # 5. AUX (contiguity metric; lower is better)
            self.aux = (self.l50 * self.n50) / self.total_length if self.total_length > 0 else 0.0

            # 6. GC content (only A/T/G/C in denominator)
            canonical = gc_count + at_count
            self.gc_content = (gc_count / canonical) * 100 if canonical > 0 else 0.0

            # 7. Additional statistics
            mean_len = self.total_length / len(self.lengths) if self.lengths else 0
            median_len = float(np.median(self.lengths)) if self.lengths else 0
            largest = self.lengths[0] if self.lengths else 0
            n_ge_1kb = sum(1 for l in self.lengths if l >= 1000)
            n_ge_10kb = sum(1 for l in self.lengths if l >= 10000)
            n_ge_100kb = sum(1 for l in self.lengths if l >= 100000)
            size_ge_1kb = sum(l for l in self.lengths if l >= 1000)
            pct_ge_1kb = (size_ge_1kb / self.total_length) * 100 if self.total_length > 0 else 0.0
            ns_total = sum(v for k, v in ambiguous_counts.items() if k == 'N')
            ns_per_100kb = (ns_total / self.total_length) * 100000 if self.total_length > 0 else 0.0
            self.ambiguous_counts = ambiguous_counts

            # 8. Build report text
            report_lines = []
            report_lines.append("--- ASSEMBLY STATISTICS ---")
            if self.min_length > 0:
                report_lines.append(f"Minimum length filter: {self.min_length:,} bp")
                report_lines.append(f"Scaffolds before filter: {total_before_filter:,} ({total_size_before_filter:,} bp)")
                report_lines.append(f"Scaffolds after filter:  {len(self.lengths):,} ({self.total_length:,} bp)")
                report_lines.append("")
            report_lines.append(f"Total Scaffolds: {len(self.lengths):,}")
            report_lines.append(f"Total Assembly Size: {self.total_length:,} bp")
            report_lines.append(f"Largest Scaffold: {largest:,} bp")
            report_lines.append(f"Mean Scaffold Length: {mean_len:,.1f} bp")
            report_lines.append(f"Median Scaffold Length: {median_len:,.1f} bp")
            report_lines.append("")
            report_lines.append(f"N50: {self.n50:,} bp (L50: {self.l50})")
            report_lines.append(f"N75: {self.n75:,} bp (L75: {self.l75})")
            report_lines.append(f"N90: {self.n90:,} bp (L90: {self.l90})")
            report_lines.append(f"AUX: {self.aux:,.1f}")
            report_lines.append("")
            report_lines.append(f"GC Content: {self.gc_content:.2f}% (canonical bases only)")
            report_lines.append(f"N bases: {ns_total:,} ({ns_per_100kb:.1f} per 100 kb)")
            report_lines.append("")
            report_lines.append("Scaffolds by size class:")
            report_lines.append(f"  >= 1 kb:   {n_ge_1kb:,}")
            report_lines.append(f"  >= 10 kb:  {n_ge_10kb:,}")
            report_lines.append(f"  >= 100 kb: {n_ge_100kb:,}")
            report_lines.append(f"  % of assembly in scaffolds >= 1 kb: {pct_ge_1kb:.2f}%")
            report_lines.append("")

            if ambiguous_counts:
                report_lines.append("Ambiguous bases detected:")
                for base, count in sorted(ambiguous_counts.items(), key=lambda x: -x[1]):
                    report_lines.append(f"  {base}: {count:,}")
                report_lines.append("")

            report_lines.append("Top 5 Scaffolds (Size in bp):")
            for i, (name, size) in enumerate(zip(self.names[:5], self.lengths[:5])):
                report_lines.append(f"  {i+1}: {name} - {size:,} bp")

            report_text = "\n".join(report_lines)

            def update_text():
                self.text_output.delete(1.0, tk.END)
                self.text_output.insert(tk.END, report_text + "\n")

            self.root.after(0, update_text)

            # 9. Save TXT report
            self.root.after(0, lambda: self.status_var.set("Waiting for TXT save location..."))
            self.save_txt_report_safe(report_text)

            # 10. Save CSV report
            self.root.after(0, lambda: self.status_var.set("Waiting for CSV save location..."))
            self.save_csv_report_safe()

            # 11. Generate plot
            self.root.after(0, lambda: self.status_var.set("Generating plot..."))
            self.generate_plot_safe()

            self.root.after(0, lambda: self.status_var.set("Done."))

        except Exception as e:
            self.root.after(0, lambda: messagebox.showerror("Error", f"An error occurred during processing:\n{str(e)}"))
            self.root.after(0, lambda: self.status_var.set("Error."))
        finally:
            self.reset_buttons()

    def reset_buttons(self):
        self.root.after(0, lambda: self.btn_process.config(state=tk.NORMAL))
        self.root.after(0, lambda: self.btn_load.config(state=tk.NORMAL))

    # ------------------------------------------------------------------
    # Save TXT
    # ------------------------------------------------------------------
    def save_txt_report_safe(self, report_text):
        save_path = filedialog.asksaveasfilename(
            title="Save TXT Report",
            defaultextension=".txt",
            filetypes=[("Text Files", "*.txt")],
            initialfile="assembly_report.txt"
        )
        if not save_path:
            return

        with open(save_path, 'w') as f:
            f.write("--- GENOME ASSEMBLY REPORT ---\n\n")
            f.write(f"Source File: {os.path.basename(self.filepath)}\n")
            f.write(report_text)
            f.write("\n\n--- SCAFFOLD LIST (Descending Order) ---\n")
            for name, size in zip(self.names, self.lengths):
                pct = (size / self.total_length) * 100 if self.total_length > 0 else 0.0
                f.write(f"{name}\t{size:,} bp\t({pct:.4f}%)\n")

        self.root.after(0, lambda: self.text_output.insert(tk.END, f"\nTXT Report saved to: {save_path}\n"))

    # ------------------------------------------------------------------
    # Save CSV
    # ------------------------------------------------------------------
    def save_csv_report_safe(self):
        save_path = filedialog.asksaveasfilename(
            title="Save CSV Report",
            defaultextension=".csv",
            filetypes=[("CSV Files", "*.csv")],
            initialfile="scaffolds.csv"
        )
        if not save_path:
            return

        with open(save_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["rank", "name", "length_bp", "percent_total"])
            for i, (name, size) in enumerate(zip(self.names, self.lengths), start=1):
                pct = (size / self.total_length) * 100 if self.total_length > 0 else 0.0
                writer.writerow([i, name, size, f"{pct:.4f}"])

        self.root.after(0, lambda: self.text_output.insert(tk.END, f"CSV Report saved to: {save_path}\n"))

    # ------------------------------------------------------------------
    # Plot (still limited to top 50 scaffolds)
    # ------------------------------------------------------------------
    def generate_plot_safe(self):
        max_scaffolds = min(len(self.lengths), 50)
        plot_lengths = self.lengths[:max_scaffolds]

        x = np.arange(1, max_scaffolds + 1)
        cumsum = np.cumsum(plot_lengths)
        green_values = np.insert(cumsum[:-1], 0, 0)
        red_values = plot_lengths

        fig, ax = plt.subplots(figsize=(10, 5), dpi=300)

        ax.bar(x, green_values, color='#2ecc71', label='Sum of Previous Scaffolds',
               edgecolor='#27ae60', linewidth=0.5)
        ax.bar(x, red_values, bottom=green_values, color='#e74c3c', label='Current Scaffold',
               edgecolor='#c0392b', linewidth=0.5)

        ax.set_xlabel('Scaffold Rank (1 = Largest)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Cumulative Size (bp)', fontsize=12, fontweight='bold')
        ax.set_title('Cumulative Bar Chart of Assembly', fontsize=14, fontweight='bold', pad=15)

        ax.grid(axis='y', linestyle='--', alpha=0.6)
        ax.set_axisbelow(True)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.legend(loc='upper left', frameon=True, facecolor='white', framealpha=0.9)

        plt.tight_layout()

        save_path = filedialog.asksaveasfilename(
            title="Save Chart",
            defaultextension=".png",
            filetypes=[("PNG Image", "*.png"), ("PDF (Vector)", "*.pdf"), ("SVG (Vector)", "*.svg")],
            initialfile="assembly_chart.png"
        )

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            self.root.after(0, lambda: self.text_output.insert(tk.END, f"Chart saved to: {save_path}\n"))

        # Display chart in the interface
        for widget in self.plot_frame.winfo_children():
            widget.destroy()

        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)


if __name__ == "__main__":
    root = tk.Tk()
    app = GenomeAnalyzerApp(root)
    root.mainloop()
