import argparse
import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import precision_recall_curve, auc
import matplotlib.pyplot as plt

# Generate k-mers from a given sequence
def generate_kmers(sequence, k=3):
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

# Load and process the dataset
def load_and_process_data(filepath, k=3):
    df = pd.read_csv(filepath, sep="\t")
    df = df[['Mature_mirna_transcript', 'Positive_Negative']]
    df['kmers'] = df['Mature_mirna_transcript'].apply(lambda x: ' '.join(generate_kmers(x, k=k)))
    return df

# Build and train the model
def train_random_forest(X_train, y_train):
    model = RandomForestClassifier(n_estimators=2, random_state=42)
    model.fit(X_train, y_train)
    return model

# Evaluate the model and plot AUC-PR
def evaluate_model(model, X_test, y_test):
    y_probs = model.predict_proba(X_test)[:, 1]  # Get probabilities for the positive class
    precision, recall, _ = precision_recall_curve(y_test, y_probs)
    aucpr = auc(recall, precision)
    print(f"AUCPR: {aucpr}")

    # Plot Precision-Recall Curve
    plt.figure()
    plt.plot(recall, precision, label=f"AUCPR = {aucpr:.2f}")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.savefig("PR_curve_random_forest.png")

# Main function
def main():
    parser = argparse.ArgumentParser(description="Train and evaluate a Random Forest model with k-mers.")
    parser.add_argument("--train", required=True, help="Path to the training dataset file.")
    parser.add_argument("--test", required=True, help="Path to the testing dataset file.")
    args = parser.parse_args()

    # Load and process data
    train_df = load_and_process_data(args.train, k=3)
    test_df = load_and_process_data(args.test, k=3)

    # Vectorize k-mers
    vectorizer = CountVectorizer()
    X_train = vectorizer.fit_transform(train_df['kmers'])
    y_train = train_df['Positive_Negative']

    X_test = vectorizer.transform(test_df['kmers'])
    y_test = test_df['Positive_Negative']

    # Train model
    model = train_random_forest(X_train, y_train)

    # Evaluate model
    evaluate_model(model, X_test, y_test)

if __name__ == "__main__":
    main()
