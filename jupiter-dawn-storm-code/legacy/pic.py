import numpy as np
import matplotlib.pyplot as plt
# Q极光数据文件路径（两列：位置(Rj)  能量(nT^2/Hz)）
q_file = r"D:\summer research\\pic\\PSD_aurora_0.9\\Q aurora.txt"  # 请修改为实际路径

# X极光数据文件路径
x_file = r"D:\summer research\\e_flux\\I-family-PSD-summary.txt"   # 请修改为实际路径

# 输出图片保存路径（可选，若不需要保存可注释掉）
output_image = r"D:\summer research\\pic\\aurora comparision\\aurora_e_I.png"
# Plotting style
q_color = 'blue'       # Color for Q aurora
x_color = 'orange'        # Color for X aurora
marker_size = 10       # Marker size
alpha = 0.7            # Transparency
# ============================================================

def read_data(file_path):
    """Read two-column data file, return position and energy arrays."""
    try:
        data = np.loadtxt(file_path)
        if data.ndim == 1:
            pos = np.array([data[0]])
            energy = np.array([data[1]])
        else:
            pos = data[:, 0]
            energy = data[:, 1]
        return pos, energy
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None, None

def main():
    # Read data
    q_pos, q_energy = read_data(q_file)
    x_pos, x_energy = read_data(x_file)
    
    if q_pos is None or x_pos is None:
        print("Data loading failed. Check file paths.")
        return
    
    # Create figure
    plt.figure(figsize=(10, 6))
    
    # Scatter plots
    plt.scatter(q_pos, q_energy, color=q_color, s=marker_size, alpha=alpha, label='Q-family aurora')
    plt.scatter(x_pos, x_energy, color=x_color, s=marker_size, alpha=alpha, label='I-family aurora')
    
    # Labels and title
    plt.xlabel(r"Radial Distance ($R_j$)", fontsize=12)
    plt.ylabel("Magnetic energy (nT²/Hz)", fontsize=12)
    plt.title("Comparison of Q and I family aurora", fontsize=14)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    
    # Save figure if output path is provided
    if output_image:
        plt.savefig(output_image, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {output_image}")
    
    plt.show()

if __name__ == "__main__":
    main()