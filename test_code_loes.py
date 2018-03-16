import lib_Loes as ll
import matplotlib.pyplot as plt

folder = r'C:\Users\Loes\Documents\GitHub\Projectstage\audiovoorbeelden\heartbeat-sounds\set_a'
data, X_list = ll.preprocess_data(folder, plot=False)

n_training_samples = data.shape[0]
n_dim = data.shape[1]

plt.figure()
plt.plot(data[:,0], data[:,1], "bx")
plt.show()

# Data plotten lukt niet, want de data heeft meer dan 2 dimenties.
# Kijk naar Xft en Xft_list
# Onset_detect doet het niet met Xft.
# Daarna verder anomaly detection uitvoeren mbv http://aqibsaeed.github.io/2016-07-17-anomaly-detection/