"""Puffin architecture (published).

PROVENANCE
----------
``Puffin`` below is the **published** architecture, transcribed from the lab's
reference implementation (``CTCF/puffin.py``). It has three sequence branches:

    conv       (4 -> 10, width 51)   the ten interpretable motif filters
    conv_inr   (4 -> 10, width 15)   the initiator branch
    conv_sim   (4 -> 32, width 3)    a trinucleotide / composition branch

Each branch is run on the sequence and on its reverse complement, softplus
activated, and deconvolved (width 601) into ten output tracks.

The only deviation from the reference file is that ``FFTConv1d`` from
``torch_fftconv`` is replaced by ``nn.Conv1d``. These are drop-in equivalent
(same parameter names and shapes, same cross-correlation semantics); checked
numerically on this cluster at max abs difference 9.5e-07 on a (1, 20, 4650)
input, against signal std 0.56. Using ``nn.Conv1d`` removes the only dependency
that was not already present in the inference environment.

This class does not load weights; see ``adapters/real_puffin.py``.
"""

from __future__ import annotations

import torch
from torch import nn

# ---------------------------------------------------------------------------
# Puffin (published)
# ---------------------------------------------------------------------------

# The ten motif filters, in the order they occupy channels 0..9 of ``conv``.
# Taken from ``motifnames_original`` in the reference implementation, which
# lists all twenty channels (0..9 = conv on the forward pass, 10..19 = conv on
# the reverse complement).
PUFFIN_MOTIF_NAMES = [
    "SP", "ETS", "CREB", "NFY", "YY1",
    "U1 snRNP", "Long Inr", "NRF1", "ZNF143", "TATA",
]

# Orientation each channel is labelled with in the reference implementation.
# Channels 0..9 carry these signs; channels 10..19 carry the opposite.
PUFFIN_MOTIF_STRAND = ["-", "-", "-", "-", "-", "+", "+", "+", "-", "-"]

# Output track layout. Puffin emits ten tracks: five assays on each strand.
PUFFIN_TRACKS_PLUS = {
    "FANTOM_CAGE": 0, "ENCODE_CAGE": 1, "ENCODE_RAMPAGE": 2,
    "GRO_CAP": 3, "PRO_CAP": 4,
}
PUFFIN_TRACKS_MINUS = {
    "FANTOM_CAGE": 9, "ENCODE_CAGE": 8, "ENCODE_RAMPAGE": 7,
    "GRO_CAP": 6, "PRO_CAP": 5,
}

# The reference ``predict`` trims this many bases from each end, because the
# 601-wide deconvolutions need that much context to be valid.
PUFFIN_TRIM = 325


class Puffin(nn.Module):
    """Published Puffin: ten motif filters, an Inr branch, and a composition branch."""

    def __init__(self):
        super().__init__()
        self.conv = nn.Conv1d(4, 10, kernel_size=51, padding=25)
        self.conv_inr = nn.Conv1d(4, 10, kernel_size=15, padding=7)
        self.conv_sim = nn.Conv1d(4, 32, kernel_size=3, padding=1)

        self.activation = nn.Softplus()
        self.softplus = nn.Softplus()

        self.deconv = nn.Conv1d(10 * 2, 10, kernel_size=601, padding=300)
        self.deconv_sim = nn.Conv1d(64, 10, kernel_size=601, padding=300)
        self.deconv_inr = nn.ConvTranspose1d(20, 10, kernel_size=15, padding=7)

        self.scaler = nn.Parameter(torch.ones(1))
        self.scaler2 = nn.Parameter(torch.ones(1))

    def forward(self, x):
        y = torch.cat([self.conv(x), self.conv(x.flip([1, 2])).flip([2])], 1)
        y_sim = torch.cat(
            [self.conv_sim(x), self.conv_sim(x.flip([1, 2])).flip([2])], 1
        )
        y_inr = torch.cat(
            [self.conv_inr(x), self.conv_inr(x.flip([1, 2])).flip([2])], 1
        )

        yact = self.activation(y)
        y_sim_act = self.activation(y_sim)
        y_inr_act = self.activation(y_inr)

        y_pred = self.softplus(
            self.deconv(yact) + self.deconv_inr(y_inr_act) + self.deconv_sim(y_sim_act)
        )
        return y_pred

    def motif_activations(self, x):
        """Post-activation motif-branch output, (batch, 20, L).

        Channels 0..9 are ``PUFFIN_MOTIF_NAMES`` at ``PUFFIN_MOTIF_STRAND``;
        channels 10..19 are the same motifs in the opposite orientation.
        Exposed so the app can show which motifs the model actually sees in a
        construct, rather than relying on a separate PWM scan.
        """
        y = torch.cat([self.conv(x), self.conv(x.flip([1, 2])).flip([2])], 1)
        return self.activation(y)
