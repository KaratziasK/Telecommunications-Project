import numpy as np
import commlib as cl

# Υπολογισμός τριγωνικού παλμού με διάρκεια TS:
# Αtri(t/T)
# p(t) = max(0, 1 - |t / (TS / 2)|), όπου p(t) είναι μηδενικό εκτός του [-TS/2, TS/2]
# t: χρόνος, TS: διάρκεια
def triangular_pulse(t, TS):
    return np.maximum(1 - np.abs(t / (TS / 2)), 0)

# Παράμετροι
p = 0 # AM: 2022120
input_bits='110101100011100010010011'
M = 2 ** (p + 3) # υπολογισμός του Μ αφού p<=5

# Έλεγχος για την είσοδο του TS = Διάρκεια συμβόλου (1 Gbit/s, δηλαδή Ts = 1 / 1GHz)
while True:
    try:
        TS = float(input("Εισάγετε διάρκεια (TS): "))
        if TS <= 0:
            print("Error: TS πρέπει να είναι θετικός")
        else:
            break
    except ValueError:
        print("Error: Εισάγετε έγκυρο αριθμό (π.χ. 1e-9).")

# Έλεγχος για την είσοδο του N = Αριθμός δειγμάτων ανά σύμβολο
while True:
    try:
        N = int(input("Εισάγετε αριθμο δειγμάτων ανά σύμβολο (N): "))
        if N <= 0:
            print("Error: N πρέπει να είναι θετικός ακαίρεος αριθμός.")
        else:
            break
    except ValueError:
        print("Error: Εισάγετε έγκυρο θετικό ακαίρεο αριθμό.")

Am = cl.pam_constellation(M)  # Δημιουργία του PAM αστερισμού για M σύμβολα
Am.set_gray_bits(m=int(np.log2(M)))  # Ορισμός Gray mapping, m = log2(M), αριθμός bit ανά σύμβολο

bits_array = np.array([int(b) for b in input_bits])  # Μετατροπή του string των bits σε array τύπου int
encoded_symbols, grouped_bits = Am.bits_to_symbols(bits_array, return_groups=True)  # Κωδικοποίηση των bits σε σύμβολα χρησιμοποιώντας τον PAM αστερισμό

duration = len(encoded_symbols) * TS  # Διάρκεια του σήματος = πλήθος συμβόλων * διάρκεια συμβόλου (8*1e-9 = 8e-9 = 8ns)
t = np.linspace(0, duration, len(encoded_symbols) * N) # Δημιουργεί τον χρονικό άξονα t (από το 0 μέχρι το duration), έχει μέσα χρονικά δείγματα
xt = np.zeros_like(t)  # Αρχικοποίηση του x(t) ως μηδενικός πίνακας ίδιου μεγέθους με τον άξονα χρόνου

# Υπολογισμός της κυματομορφής x(t) = Σk akp(t - kTs))
for k, a in enumerate(encoded_symbols):  
    # Για κάθε σύμβολο του PAM αστερισμού:
    pulse = triangular_pulse(t - k * TS, TS) # Υπολογισμός τριγωνικού παλμού p(t - kTs) για τη θέση του τρέχοντος συμβόλου
    xt += a * pulse # Πρόσθεση του βάρους του παλμού με βάση το σύμβολο a (x(t) = Σk akp(t - kTs))

# Γραφική παράσταση του x(t)
cl.plot_signal(t * 1e9, xt, plot_type='-', close_all=False,  
               xlabel="t", ylabel="x(t)",  
               title="Κυματομορφη PAM με τριγωνικους παλμους", 
               show_grid=True)

# Εκτύπωση αποτελεσμάτων
print("Εισαγόμενα Bits:", input_bits)  
# Αρχική ακολουθία bits
print("Ομαδοποιημένα Bits:", grouped_bits)  
# Bits ομαδοποιημένα σε ομάδες μεγέθους log2(M)
print("Κωδικοποιημένα Σύμβολα:", encoded_symbols)  
# Τα σύμβολα που αντιστοιχούν στην ακολουθία bits



#============================================================
# Part 2

sigma2_a = 1.0  # διακύμανση των συμβόλων
t = cl.time_axis(0, len(encoded_symbols) * TS, len(encoded_symbols) * N)

# Δημιουργία σήματος με το x(t) και τον χρονικό άξονα
signal = cl.signal(t=t, samples=xt) # Υπολογισμός Fourier Φάσματος για την αριθμητική φασματική πυκνότητα
signal.calc_spectrum() # Υπολογισμός Fourier Φάσματος για την αριθμητική φασματική πυκνότητα  {X(f) = ∫ x(t) e^(-j2πft) dt}

# Υπολογισμός της φασματικής πυκνότητας ισχύος από το φάσμα
SX_numerical = signal.power_density() # Υπολογίζει: SX(f) = (1/T) * |X(f)|^2, όπου T είναι η διάρκεια του σήματος
f = signal.set_frequency_axis() # Δημιουργεί άξονα συχνοτήτων f = [-Fs/2, Fs/2), όπου Fs = 1/Δt, με Δt το διάστημα δειγματοληψίας


# Μαθηματικός υπολογισμός της φασματικής πυκνότητας ισχύος
f_math = np.linspace(-1 / (2 * TS), 1 / (2 * TS), len(SX_numerical)) # Δημιουργεί άξονα συχνοτήτων [-1/(2TS), 1/(2TS)]
SX_math = (sigma2_a / TS) * (np.sinc(f_math * TS))**2 # Υπολογίζει: SX(f) = (σ^2_a / TS) * sinc^2(f * TS), όπου sinc(x) = sin(πx) / (πx)

# Γραφική Αναπαράσταση
cl.plot_signal(f * 1e-9, SX_numerical, 
               xlabel='f[GHz]', 
               ylabel='SX', 
               title='Αριθμητική Εκτίμηση Φασματικής Πυκνότητας Ισχύος', 
               show_grid=True, close_all=False)

cl.plot_signal(f_math * 1e-9, SX_math, 
               xlabel='f*TS', 
               ylabel='SX', 
               title='Μαθηματικός Υπολογισμός Φασματικής Πυκνότητας Ισχύος', 
               show_grid=True,close_all=False)
