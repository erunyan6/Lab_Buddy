package sanity_check

import (
	"fmt"

	"lab_buddy_go/config"
)

// Simple sanity check
func Run(args []string) {
	fmt.Printf("Successfully running Lab Buddy! (%s)\n", version_control.Main_version)
}