// benchmark.go
// A reusable benchmarking module for Lab Buddy
// Measures execution time and memory usage for any wrapped function

package benchmark

import (
	"fmt"
	"os"
	"runtime"
	"time"
)

// Run wraps any function to measure its runtime and memory usage.
func Run(label string, f func()) {
	fmt.Printf("[Benchmark] Running: %s\n", label)

	// Snapshot environment info
	fmt.Println("[Benchmark] Timestamp:", time.Now().Format(time.RFC1123))
	host, err := os.Hostname()
	if err == nil {
		fmt.Println("[Benchmark] Hostname:", host)
	}
	fmt.Println("[Benchmark] Go Version:", runtime.Version())
	fmt.Printf("[Benchmark] OS/Arch: %s/%s\n", runtime.GOOS, runtime.GOARCH)

	// Prepare for benchmark
	runtime.GC()
	var memStart, memEnd runtime.MemStats
	runtime.ReadMemStats(&memStart)
	start := time.Now()
	numCPU := runtime.NumCPU()
	startGoroutines := runtime.NumGoroutine()

	// Run benchmarked function
	f()

	elapsed := time.Since(start)
	runtime.ReadMemStats(&memEnd)
	endGoroutines := runtime.NumGoroutine()

	// Report resource usage
	fmt.Printf("[Benchmark] Time Elapsed: %v\n", elapsed)
	fmt.Printf("[Benchmark] Memory Used: %.2f MB\n", float64(memEnd.Alloc-memStart.Alloc)/1024.0/1024.0)
	fmt.Printf("[Benchmark] Total Allocated: %.2f MB\n", float64(memEnd.TotalAlloc-memStart.TotalAlloc)/1024.0/1024.0)
	fmt.Printf("[Benchmark] Peak Heap: %.2f MB\n", float64(memEnd.HeapAlloc)/1024.0/1024.0)
	fmt.Printf("[Benchmark] GC Cycles: %d\n", memEnd.NumGC-memStart.NumGC)
	fmt.Printf("[Benchmark] CPU Cores: %d\n", numCPU)
	fmt.Printf("[Benchmark] Goroutines Started: %d → %d\n", startGoroutines, endGoroutines)
	fmt.Println("[Benchmark] ----------------------------------------")
}
