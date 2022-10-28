import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from math import sqrt

fineza_coef_alph_bet = 50

partition = np.arange(0, fineza_coef_alph_bet)
alphas = np.round(np.sqrt(partition/(partition.shape[0]-1)), 4)
betas  = np.round(np.sqrt(1 - np.power(alphas, 2)), 4)

fig = plt.figure(figsize=(8, 6))
gs = gridspec.GridSpec(nrows = 2,ncols=1, figure=fig, height_ratios=(5, 5),
                left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.05, hspace=0)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0])

alpha_2, beta_2 = [alphas[i]**2 for i in range(len(alphas))], [betas[i]**2 for i in range(len(betas))]
ax1.bar(partition, alpha_2, width=0.5, color='red')
ax1.bar(partition, beta_2, bottom=alpha_2, width=0.5)
ax2.bar(partition, alphas, width=0.5, color='orange')
ax2.bar(partition, betas, bottom=alphas, width=0.5)
ax2.hlines(1,xmin=0,xmax=50, color= 'dark')
plt.show()
