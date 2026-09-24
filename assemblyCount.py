import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import os

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
        self.root.geometry("1000x750")
        self.root.configure(bg="#ecf0f1")

        # ttk Styling (Blue tones)
        self.style = ttk.Style()
        self.style.theme_use('clam')
        
        # Base colors
        self.bg_color = "#ecf0f1"  # Light blue-gray
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

        self.create_widgets()
        self.filepath = None
        self.lengths = []
        self.total_length = 0
        self.gc_content = 0.0
        self.n50, self.l50 = 0, 0
        self.n90, self.l90 = 0, 0

    def create_widgets(self):
        # Top Frame (Buttons and Status)
        top_frame = ttk.Frame(self.root, padding="10")
        top_frame.pack(fill=tk.X)

        ttk.Label(top_frame, text="assemblyCount", style='Header.TLabel').pack(side=tk.LEFT, padx=10)

        self.btn_load = ttk.Button(top_frame, text="Load FASTA", command=self.load_fasta)
        self.btn_load.pack(side=tk.RIGHT, padx=5)

        self.btn_process = ttk.Button(top_frame, text="Process & Save Results", command=self.process_data, state=tk.DISABLED)
        self.btn_process.pack(side=tk.RIGHT, padx=5)

        # Middle Frame (Text Information)
        middle_frame = ttk.Frame(self.root, padding="10")
        middle_frame.pack(fill=tk.BOTH, expand=False)

        info_frame = ttk.LabelFrame(middle_frame, text=" Assembly Results ", padding="10")
        info_frame.pack(fill=tk.X)

        self.text_output = tk.Text(info_frame, height=10, width=80, bg="white", fg=self.text_color, font=('Courier', 10), relief="flat", borderwidth=1)
        self.text_output.pack(fill=tk.X)

        # Bottom Frame (Plot)
        self.plot_frame = ttk.LabelFrame(self.root, text=" Cumulative Bar Chart ", padding="10")
        self.plot_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

    def load_fasta(self):
        self.filepath = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA Files", "*.fasta *.fa *.fna"), ("All Files", "*.*")]
        )
        if self.filepath:
            self.btn_process.config(state=tk.NORMAL)
            self.text_output.delete(1.0, tk.END)
            self.text_output.insert(tk.END, f"File loaded: {os.path.basename(self.filepath)}\n")
            self.text_output.insert(tk.END, "Click 'Process & Save Results' to start the analysis.\n")

    def parse_fasta(self):
        """Reads the FASTA file and returns scaffold lengths and GC count."""
        lengths = []
        gc_count = 0
        total_bases = 0
        current_len = 0
        
        with open(self.filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_len > 0:
                        lengths.append(current_len)
                        current_len = 0
                else:
                    current_len += len(line)
                    total_bases += len(line)
                    gc_count += line.upper().count('G') + line.upper().count('C')
            if current_len > 0:
                lengths.append(current_len)
                
        # Sort in descending order of length
        lengths.sort(reverse=True)
        return lengths, total_bases, gc_count

    def calculate_nx_lx(self, lengths, total_len, x):
        """Calculates Nx and Lx (e.g., N50 and L50)."""
        target = total_len * (x / 100.0)
        cum_sum = 0
        for i, length in enumerate(lengths):
            cum_sum += length
            if cum_sum >= target:
                return length, (i + 1)
        return 0, 0

    def process_data(self):
        if not self.filepath:
            return

        try:
            self.text_output.insert(tk.END, "Processing data...\n")
            self.root.update()

            # 1. Read FASTA
            self.lengths, total_bases, gc_count = self.parse_fasta()
            
            if not self.lengths:
                messagebox.showerror("Error", "No sequences found in the FASTA file.")
                return

            # 2. Total assembly size
            self.total_length = sum(self.lengths)

            # 3 & 4. N50, L50, N90, L90
            self.n50, self.l50 = self.calculate_nx_lx(self.lengths, self.total_length, 50)
            self.n90, self.l90 = self.calculate_nx_lx(self.lengths, self.total_length, 90)

            # 5. C/G content
            self.gc_content = (gc_count / total_bases) * 100 if total_bases > 0 else 0.0

            # Display on screen
            self.text_output.delete(1.0, tk.END)
            self.text_output.insert(tk.END, f"--- ASSEMBLY STATISTICS ---\n")
            self.text_output.insert(tk.END, f"Total Scaffolds: {len(self.lengths)}\n")
            self.text_output.insert(tk.END, f"Total Assembly Size: {self.total_length:,} bp\n")
            self.text_output.insert(tk.END, f"N50: {self.n50:,} bp (L50: {self.l50})\n")
            self.text_output.insert(tk.END, f"N90: {self.n90:,} bp (L90: {self.l90})\n")
            self.text_output.insert(tk.END, f"C/G Content: {self.gc_content:.2f}%\n")
            self.text_output.insert(tk.END, f"\nTop 5 Scaffolds (Size in bp):\n")
            for i, size in enumerate(self.lengths[:5]):
                self.text_output.insert(tk.END, f"  {i+1}: {size:,}\n")

            # Save TXT file
            self.save_txt_report()

            # Generate Plot
            self.generate_plot()

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during processing:\n{str(e)}")

    def save_txt_report(self):
        """Saves information 1 to 5 (including L50/L90) to a TXT file."""
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
            f.write(f"Total Scaffolds: {len(self.lengths)}\n")
            f.write(f"Total Assembly Size: {self.total_length:,} bp\n")
            f.write(f"N50: {self.n50:,} bp\n")
            f.write(f"L50: {self.l50}\n")
            f.write(f"N90: {self.n90:,} bp\n")
            f.write(f"L90: {self.l90}\n")
            f.write(f"C/G Content: {self.gc_content:.2f}%\n\n")
            f.write("--- SCAFFOLD LIST (Descending Order) ---\n")
            for i, size in enumerate(self.lengths):
                f.write(f"Scaffold {i+1}: {size:,} bp\n")
        
        self.text_output.insert(tk.END, f"\nTXT Report saved to: {save_path}\n")

    def generate_plot(self):
        """Generates the cumulative bar chart as specified."""
        # Limit the number of scaffolds plotted to avoid clutter (e.g., Top 50)
        max_scaffolds = min(len(self.lengths), 50)
        plot_lengths = self.lengths[:max_scaffolds]
        
        # Cumulative calculation
        x = np.arange(1, max_scaffolds + 1)
        cumsum = np.cumsum(plot_lengths)
        
        # Green bars: Sum of all previous scaffolds (for the first bar, it is 0)
        green_values = np.insert(cumsum[:-1], 0, 0)
        
        # Red bars: The current scaffold
        red_values = plot_lengths

        # Create high-quality figure
        fig, ax = plt.subplots(figsize=(10, 5), dpi=300)
        
        # Plot bars
        ax.bar(x, green_values, color='#2ecc71', label='Sum of Previous Scaffolds', edgecolor='#27ae60', linewidth=0.5)
        ax.bar(x, red_values, bottom=green_values, color='#e74c3c', label='Current Scaffold', edgecolor='#c0392b', linewidth=0.5)

        # Aesthetic configurations for the plot
        ax.set_xlabel('Scaffold Rank (1 = Largest)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Cumulative Size (bp)', fontsize=12, fontweight='bold')
        ax.set_title('Cumulative Bar Chart of Assembly', fontsize=14, fontweight='bold', pad=15)
        
        # Add soft grid on Y axis
        ax.grid(axis='y', linestyle='--', alpha=0.6)
        ax.set_axisbelow(True)
        
        # Remove top and right borders
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Add legend
        ax.legend(loc='upper left', frameon=True, facecolor='white', framealpha=0.9)
        
        plt.tight_layout()

        # Save figure in high quality
        save_path = filedialog.asksaveasfilename(
            title="Save Chart",
            defaultextension=".png",
            filetypes=[("PNG Image", "*.png"), ("PDF (Vector)", "*.pdf"), ("SVG (Vector)", "*.svg")],
            initialfile="assembly_chart.png"
        )
        
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            self.text_output.insert(tk.END, f"Chart saved to: {save_path}\n")

        # Display chart on the interface
        for widget in self.plot_frame.winfo_children():
            widget.destroy()
            
        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

if __name__ == "__main__":
    root = tk.Tk()
    app = GenomeAnalyzerApp(root)
    root.mainloop()
