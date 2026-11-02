#!/usr/bin/env bash
set -eu
meter="${1:-$(pwd)/build/src/hycu-benchmark-measure}"
fixture=$(mktemp -d)
"$meter" --timeout-seconds 0 "$fixture/unlimited.time" "$fixture/unlimited.log" /bin/sleep 2
grep -q 'exit_code=0' "$fixture/unlimited.time"
"$meter" "$fixture/default.time" "$fixture/default.log" /bin/true
grep -q 'exit_code=0' "$fixture/default.time"
if "$meter" --timeout-seconds 1 "$fixture/timeout.time" "$fixture/timeout.log" /bin/sleep 10; then
    exit 1
else
    test "$?" -eq 124
fi
grep -q 'exit_code=124' "$fixture/timeout.time"
if "$meter" --timeout-seconds -1 "$fixture/invalid.time" "$fixture/invalid.log" /bin/true; then exit 1; fi
test ! -e "$fixture/invalid.time"
# Existing outputs must not be overwritten.
if "$meter" "$fixture/default.time" "$fixture/overwrite.log" /bin/true; then exit 1; fi
test ! -e "$fixture/overwrite.log"
# Terminating the meter must terminate its detached child group too.
"$meter" --timeout-seconds 0 "$fixture/signal.time" "$fixture/signal.log" /bin/sh -c 'echo $$; exec /bin/sleep 30' &
meter_pid=$!
for i in $(seq 1 100); do
    test -s "$fixture/signal.log" && break
    sleep 0.02
done
child_pid=$(head -n 1 "$fixture/signal.log")
test -n "$child_pid"
kill -TERM "$meter_pid"
if wait "$meter_pid"; then exit 1; else test "$?" -eq 143; fi
if kill -0 "$child_pid" 2>/dev/null; then exit 1; fi
grep -q 'exit_code=143' "$fixture/signal.time"
echo "PASS: default/unlimited/finite timeouts, invalid argument, overwrite protection, signal cleanup. Fixtures: $fixture"
