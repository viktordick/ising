import os

Command(
    'target/release/analyze',
    ['src/bin/analyze.rs', 'Cargo.toml'],
    'cargo build --release --bin analyze',
)

nc = Environment()
nc.Decider('MD5-timestamp')
results = Glob('result/*/*', strings=True)
if os.path.exists('data'):
    for datafile in Glob('data/*/*', strings=True):
        results.append(
            nc.Command(
                'result'+datafile[4:], 
                [datafile, 'target/release/analyze'], 
                'target/release/analyze 100 0.1 $SOURCE > $TARGET',
            )
        )
